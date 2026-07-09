function [biochemFluid, model, schedule, state0] = setupH2StorageExample(varargin)
% Set up a detailed H2 storage case for bio‑reactive transport studies
% with a cyclic injection/production schedule.
%
% SYNOPSIS:
%   [biochemFluid, model, schedule, state0] = setupH2StorageExample()
%   [biochemFluid, model, schedule, state0] = setupH2StorageExample('bacteriamodel', true, ...)
%
% OPTIONAL PARAMETERS (property/value pairs):
%   bacteriamodel       - Enable bacterial growth/decay (default: false)
%   bactDiffusion       - Enable microbial diffusion (default: false)
%   chemotaxisEffect    - Enable bacterial chemotaxis (default: false)
%   molecularDiffusion  - Enable molecular diffusion (default: false)
%   molecularDispersion - Enable mechanical dispersion (default: false)
%   bioClogging         - Enable bio‑clogging (default: false)
%   nbact0              - Initial bacterial concentration (default: 1e-6)
%   ncycles             - Number of injection/production cycles (default: 1)
%
% RETURNS:
%   biochemFluid - Biochemical reaction definition
%   model        - BiochemistryModel
%   schedule     - Cyclic schedule (build‑up, rest, injection, idle, production, idle)
%   state0       - Initial state

    require ad-props compositional deckformat h2-biochem

    % --- Parse optional parameters ---
    opt = struct('bacteriamodel', false, ...
                 'bactDiffusion', false, ...
                 'chemotaxisEffect', false, ...
                 'molecularDiffusion', false, ...
                 'molecularDispersion', false, ...
                 'bioClogging', false, ...
                 'nbact0', 10, ...
                 'ncycles', 1);
    opt = merge_options(opt, varargin{:});

    %% 1. Grid and rock (use the simplecomp deck for a real 3D grid)
    pth = getDatasetPath('simplecomp');
    fn  = fullfile(pth, 'SIMPLE_COMP.DATA');
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    G = initEclipseGrid(deck);
    G = computeGeometry(G);
    rock = initEclipseRock(deck);
    rock = compressRock(rock, G.cells.indexMap);

    %% 2. Fluid components and EOS
    % Five components: H2O, H2, CO2, CH4, CH3COOH
    compFluid = TableCompositionalMixture( ...
        {'Water', 'Hydrogen', 'CarbonDioxide', 'Methane', 'AceticAcid'}, ...
        {'H2O', 'H2', 'CO2', 'C1', 'CH3COOH'});
    biochemFluid = TableBioChemMixture({'MethanogenicArchae', 'AcetogenicBacteria'}, ...
                                       {'bactM', 'bactA'});
    eos = SoreideWhitsonEos(G, compFluid, 'msalt', 0);  % fresh water

    % Simple fluid (used as a base for transport properties)
    fluid = initSimpleADIFluid('phases', 'OG', ...
        'mu',  [1.3059*centi*poise, 0.01763*centi*poise], ...
        'rho', [999.7, 1.2243] .* kilogram/meter^3, ...
        'pRef', 100*barsa, ...
        'c',   [5.0e-5/barsa, 1.0/barsa], ...
        'n',   [2, 2], ...
        'smin',[0.2, 0.05]);
    Pe = 0.1*barsa;
    fluid.pcOG = @(sg) Pe * max((1 - sg - 0.2) ./ (1 - 0.2), 1e-5).^(-1/2);

    %% 3. Assemble model
    backend = DiagonalAutoDiffBackend('modifyOperators', true);
    model = BiochemistryModel(G, rock, fluid, compFluid, biochemFluid, ...
        true, backend, ...   % explicit = true
        'water', false, 'oil', true, 'gas', true, ...
        'bacteriamodel', opt.bacteriamodel, ...
        'bactDiffusion', opt.bactDiffusion, ...
        'chemotaxisEffect', opt.chemotaxisEffect, ...
        'molecularDiffusion', opt.molecularDiffusion, ...
        'molecularDispersion', opt.molecularDispersion, ...
        'liquidPhase', 'O', 'vaporPhase', 'G');
    model.EOSModel = eos;

    % Bio‑clogging (if requested)
    if opt.bioClogging && opt.bacteriamodel
        nc = [180,180];               % critical bacteria concentration
        cp = [0.5, 0.5];               % scaling coefficient
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        model = setupBioCloggingModel(model, nbact0, nc, cp, true);
    else
        % Attach a dummy bio‑clogging model to avoid errors
        nc = [180,180];               % critical bacteria concentration
        cp = [0.0, 0.0];               % scaling coefficient
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        model = setupBioCloggingModel(model, nbact0, nc, cp, false);
    end

    %% 4. Initial state
    p0 = 82*barsa;                     % initial pressure
    T0 = 273.15 + 40;                  % 40 °C
    s0 = [0.21, 0.79];                   % fully water‑saturated
    z0_overall = [0.8480, 0.5e-5, 1.0e-5, 0.1530,0.5e-5];  % pure water (mole fractions)
    % Replicate to all cells
    ncell = G.cells.num;
    z0 = repmat(z0_overall, ncell, 1);

    if opt.bacteriamodel
        nbact0 = opt.nbact0 * ones(1, biochemFluid.nbioreact);
        state0 = initCompositionalStateBacteria(model, p0, T0, s0, z0, nbact0, ...
                                                 eos);
    else
        state0 = initCompositionalState(G, p0, T0, s0, z0, eos);
    end

    %% 5. Wells (CO2‑rich injector, H2‑rich injector, producer)
    injCell = 1;
    prodCell = G.cells.num;
    W = addWell([], G, rock, injCell, ...
                'Type', 'rate', ...        % controlled by surface rate
                'Val', 1 * kilogram/day, ... % injection rate (reservoir volume)
                'components',  [0.0, 0.95, 0.05, 0.0, 0.0],...
                'Radius', 0.1, ...
                'Comp_i', [0, 1], ...      % inject pure H2 (mole frac)
                'Name', 'Injector');
    W = addWell(W, G, rock, prodCell, ...
                'Type', 'bhp', ...         % controlled by bottom-hole pressure
                'Val', 20*barsa, ...
                'Radius', 0.1, ...
                'Comp_i', [0, 1], ...      % inject pure H2 (mole frac)
                'components',  [0.0, 0.95, 0.05, 0.0, 0.0],...
                'Name', 'Producer');

    %% 7. Create a three-stage schedule: inject -> shut-in -> produce
    % Stage lengths (days)
    injDays  = 30*day;
    shutDays = 30*day;
    prodDays = 30*day;

    % Timestep size (all stages use 1-day steps)
    dt = 1*day;

    % Control definitions
    % Control 1: injection (both wells active)
    ctrl(1).W = W;                    % as defined (injector on, producer on)
    % Control 2: shut-in (both wells inactive)
    ctrl(2).W = W;
    ctrl(2).W(1).status = false;      % shut injector
    ctrl(2).W(2).status = false;      % shut producer
    % Control 3: production (injector off, producer on)
    ctrl(3).W = W;
    ctrl(3).W(1).status = false;      % injector shut
    ctrl(3).W(2).status = true;       % producer open (same BHP as before)

    % Build step arrays
    nStepsInj  = ceil(injDays / day);
    nStepsShut = ceil(shutDays / day);
    nStepsProd = ceil(prodDays / day);
    nTotal     = nStepsInj + nStepsShut + nStepsProd;

    step.val     = repmat(dt, nTotal, 1);
    step.control = [ones(nStepsInj, 1); ...
                    2*ones(nStepsShut, 1); ...
                    3*ones(nStepsProd, 1)];

    schedule.control = ctrl;
    schedule.step    = step;

    % Remove unnecessary fields (safety)
    for i = 1:numel(schedule.control(1).W)
        schedule.control(1).W(i).lims = [];
    end
end