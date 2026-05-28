%% 2D Compositional Hydrogen Storage with Bio‑Clogging: Diffusion & Dispersion Analysis
% ===========================================================================
% This script compares four scenarios of a bio‑clogged hydrogen storage
% simulation with bacteria effects enabled:
%   1) No molecular diffusion, no mechanical dispersion (pure advection)
%   2) Molecular diffusion only
%   3) Mechanical dispersion only
%   4) Both molecular diffusion and mechanical dispersion
%
% The aim is to isolate the impact of these two transport mechanisms on
% the evolution of hydrogen, CO2, CH4 distributions, and on biomass‑induced
% porosity/permeability changes.
% ===========================================================================

mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr test-suite spe10
mrstModule add compositional h2-biochem h2store

%% 1. Load and setup the base model (common to all scenarios)
baseName = 'H2_STORAGE_DOME_TRAP';
dataPath = getDatasetPath('h2storage');
dataFile = fullfile(dataPath, 'H2STORAGE_RS.DATA');
deck = readEclipseDeck(dataFile);

% Black‑oil model and schedule (10 injection cycles)
[~, ~, state0Bo, modelBo, scheduleBo, ~] = modelForSimple2DAquifer(deck, 'numcycles', 10);

% Convert to compositional model
model0 = convertBlackOilModelToCompositionalModel(modelBo);
state0 = convertBlackOilStateToCompositional(modelBo, state0Bo);

% Fluid and EOS
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});
EOS = SoreideWhitsonEos([], compFluid);
model0.EOSModel = EOS;

nc = model0.G.cells.num;
T0 = 273.15 + 44.35;  % Initial temperature [K]
comp0 = repmat([0.8480, 1.0e-5, 1.0e-5, 0.1530], nc, 1);

%% 2. Bio‑clogging parameters (identical for all four cases)
bacteriamodel = true;       % bacteria always present
clogModel     = true;       % clogging always active
nbact0 = 5;                 % initial bacteria concentration
nc_bact = 120;              % critical concentration
cp = 1.0;                   % clogging coefficient

% Setup base bio‑clogging model and store original porosity/permeability
[model0, poro0, perm0] = setupBioCloggingModel(model0, nbact0, nc_bact, cp);

%% 3. Update schedule controls for compositional injection
schedule = scheduleBo;
for i = 1:numel(schedule.control)
    schedule.control(i).W.compi = [0, 1];  % well components (water, oil)
    if strcmp(schedule.control(i).W.name, 'cushion') && i < 11
        schedule.control(i).W.components = [0.0, 0.1, 0.9, 0.0];
    else
        schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
    end
    schedule.control(i).W.T = T0;          % injection temperature
    schedule.control(i).bc = [];           % remove boundary conditions
end

%% 4. Common initial state (with bacteria)
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);

%% 5. Helper to create BiochemistryModel with given diffusion/dispersion flags
createModel = @(doDiff, doDisp) BiochemistryModel( ...
    model0.G, model0.rock, model0.fluid, compFluid, ...
    false, diagonal_backend, ...               % false = don't re‑initialise from deck
    'oil', true, 'gas', true, ...
    'bacteriamodel', true, ...
    'molecularDiffusion', doDiff, ...
    'molecularDispersion', doDisp, ...
    'liquidPhase', 'O', 'vaporPhase', 'G');


% 6.1 No diffusion, no dispersion
model_00 = createModel(false, false);
%% 6. Solve the four scenarios (all with bacteria & clogging)
nls = NonLinearSolver();
nls.LinearSolver = selectLinearSolverAD(model_00); % use same linear solver for all
state0 = initCompositionalStateBacteria(model_00, state0.pressure, T0, state0.s, comp0, nbact0, EOS);
prob_00  = packSimulationProblem(state0, model_00, schedule, ...
    [baseName '_noDiff_noDisp'], 'NonLinearSolver', nls);
simulatePackedProblem(prob_00);
[ws00, states00] = getPackedSimulatorOutput(prob_00);

% % 6.2 Diffusion only
model_D  = createModel(true, false);
prob_D   = packSimulationProblem(state0, model_D, schedule, ...
    [baseName '_diffOnly'], 'NonLinearSolver', nls);
simulatePackedProblem(prob_D,'RestartStep',123);
[wsD, statesD] = getPackedSimulatorOutput(prob_D);
% 
% % 6.3 Dispersion only
model_Dp = createModel(false, true);
prob_Dp  = packSimulationProblem(state0, model_Dp, schedule, ...
    [baseName '_dispOnly'], 'NonLinearSolver', nls);
simulatePackedProblem(prob_Dp);
[wsDp, statesDp] = getPackedSimulatorOutput(prob_Dp);

% 6.4 Both diffusion and dispersion
model_DD = createModel(true, true);
prob_DD  = packSimulationProblem(state0, model_DD, schedule, ...
    [baseName '_bothDiffDisp'], 'NonLinearSolver', nls);
simulatePackedProblem(prob_DD, 'RestartStep',1);
[wsDD, statesDD] = getPackedSimulatorOutput(prob_DD);

%% 7. Visualise results – states (toolbar plots)
figure; plotToolbar(model0.G, states00); title('No Diffusion, No Dispersion');
figure; plotToolbar(model0.G, statesD);  title('Diffusion Only');
figure; plotToolbar(model0.G, statesDp); title('Dispersion Only');
figure; plotToolbar(model0.G, statesDD); title('Both Diffusion & Dispersion');

%% 8. Visualise well outputs – comparison
% Collect well solutions for all scenarios (assume 2 wells: injector=1, producer=2)
% Use the same well index (e.g., 1 = injector) for plotting. You can change index.
ws_combined = {ws00, wsD, wsDp, wsDD};
figure; plotWellSols(ws_combined, cumsum(schedule.step.val)/day, ...
    'datasetnames', {'No Diff/No Disp', 'Diffusion Only', 'Dispersion Only', 'Both Diff & Disp'});

%% 9. Optional: Porosity comparison at final state
figure;
subplot(2,2,1); plotCellData(model0.G, states00{end}.rock.poro); title('No Diff/No Disp');
subplot(2,2,2); plotCellData(model0.G, statesD{end}.rock.poro);  title('Diffusion Only');
subplot(2,2,3); plotCellData(model0.G, statesDp{end}.rock.poro); title('Dispersion Only');
subplot(2,2,4); plotCellData(model0.G, statesDD{end}.rock.poro); title('Both Diff & Disp');
sgtitle('Final Porosity');

%% Copyright notice
% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% …
% </html>