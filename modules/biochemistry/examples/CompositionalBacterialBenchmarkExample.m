clear all;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% Read the Eclipse deck file containing the simulation data
% Change input fil by UHS_BENCHMARK_RS_SALT.DATA for SALT EFFECTS
%deck = readEclipseDeck('/home/elyes/Documents/mrst-2023b/spe11-utils/deck_H2/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');
deck = readEclipseDeck('/home/elyes/Documents/mrst-2023b/spe11-utils/TUC_UHS_Benchmark/Simulation Cases/UHS_Benchmark_LowH2.DATA');

%% Prepare simulation parameters and initial state
[~, options, state0, model, schedule, ~] = modified_uhs_benchmark_compositional(deck);




%%=========initialisation proprietes du fluide
%===Compositional fluid model (initialization with the CoolProp library)=
compFluid =TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide',...
    'Nitrogen','Methane'}, {'H2O','H2','CO2','N2','CH4'});

% %Brine-Gas (H2)
% [rhow,rhog]=deal(1000* kilogram/meter^3,8.1688* kilogram/meter^3); %density kilogram/meter^3;
% [viscow,viscog]=deal(1.0*centi*poise,0.0094234*centi*poise);%viscosity
% [cfw,cfg]=deal(0,8.1533e-3/barsa); %compressibility
%  
% [srw,src]=deal(0.0,0.0);
% Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% [Pmaxz,Pref1,Pminz,Pe]=deal(95*barsa,114*barsa,120*barsa,0.1*barsa); %pressions
% 
% %initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
% fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
%                          'rho',[rhow,rhog],'pRef',Pref1,...
%                          'c',[cfw,cfg],'n',[2,2],'smin',[srw,src]);
% 
% % Pression capillaire
% pcOG = @(so) Pe * so.^(-1/2);
% fluid.pcOG = @(sg) pcOG(max((1-sg-srw)./(1-srw), 1e-5)); %@@
% 

% T = 100*day;
% pv=sum(poreVolume(G,rock))/T;
% rate = 100*pv;
% niter = 30;

% %%==============Conditions aux limites et puits============
% bc = [];
% %puit d'injection
% W = [];
% W = verticalWell(W, G, rock,1,1,nz, 'comp_i', [0, 1],'Radius',0.5,...
%     'name', 'Injector', 'type', 'rate','Val',rate, 'sign', 1);
% W(1).components = [0.0 0.95 0.05 0.0 0.0];
% % for i = 1:numel(W)
% %     W(i).components = info.injection;
% % end
%%==============model compositionnal================
% arg = {model.G, model.rock, model.fluid, compFluid,...
%     'water', false, 'oil', true, 'gas', true,... % water-oil system
% 	'bacteriamodel', true,'diffusioneffect',false,'liquidPhase', 'O',...
%     'vaporPhase', 'G'}; % water=liquid, gas=vapor
% model = BioChemsitryModel(arg{:});
model.outputFluxes = false;
% %===Conditions initiales=====================
% T0=317.5;
% s0= [0.8 0.2]; %initial saturations  Sw=1
% z0 = [0.8,0.0,0.006,0.018,0.176]; %initial composition: H2O,H2,CO2,N2,CH4.
% 
% %===Bacteria model===========================
% if model.bacteriamodel
%     nbact0=10^6;
%     state0 = initCompositionalStateBacteria(model,Phydro0,T0,s0,z0,nbact0);
% else
%     state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
% end
% 
% 
% %===Ajout d'un terme source====================
% src=[];

%===Resolution pression/transport========================
% deltaT = T/niter;
% schedule = simpleSchedule(repmat(deltaT,1,niter),'bc', bc,'src', src,'W',W);
% nls = NonLinearSolver('useRelaxation', true);
%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

name = 'UHS_BENCHMARK_COMPOSITIONAL';
%% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Run the simulation
simulatePackedProblem(problem,'restartStep',1);
%% gGet reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
%===Plottings============================================
time=0;
figure;
for i= 1:niter
    x = G.cells.centroids(:,1);
    z = G.cells.centroids(:,3);
    X = reshape(x, [nx,nz]);
    Z = reshape(z, [nx,nz]);
    zH2 = reshape(states{i}.components(:,2), [nx,nz]);
    zH2O = reshape(states{i}.components(:,1), [nx,nz]);
    zCO2 = reshape(states{i}.components(:,3), [nx,nz]);
    Sw = reshape(states{i}.s(:,1), [nx,nz]);
    Sg = reshape(states{i}.s(:,2), [nx,nz]);
    Pres= reshape(states{i}.pressure, [nx,nz]);
    nbacteria=reshape(states{i}.nbact, [nx,nz]);

    subplot(2,2,1);   
    contourf(X,Z,Sw,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Water saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,2,2);   
    contourf(X,Z,Sg,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Gas saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,2,3); 
    contourf(X,Z,nbacteria,60,'EdgeColor','auto');
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('nbact','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar

    subplot(2,2,4); 
    contourf(X,Z,zH2,60,'EdgeColor','auto');
    clim([0 0.8])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('z_H2','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar
   
    time = time + deltaT;
    title(sprintf('injection duration = %.2f days',convertTo(time,day)))
    pause(0.001)
end

    function [description, options, state0, model, schedule, plotOptions] = modified_uhs_benchmark_compositional(deck,varargin)
% Example modified from:
% Hogeweg, S., Strobel, G., & Hagemann, B. (2022). Benchmark study for the simulation of underground 
% hydrogen storage operations. Comput Geosci, 26, 1367–1378. 
% https://doi.org/10.1007/s10596-022-10163-5
%
% For further details on this case, refer to:
% Ahmed, E., et al. (2024). Phase behavior and black-oil simulations of hydrogen storage in saline aquifers. 
% Advances in Water Resources, 191, 104772.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
% Step 1: Test case description and options
% Description
description = ['Modified benchmark for hydrogen storage with ', ...
               'multiple injection/production cycles using the black-oil model'];

%% Optional input arguments: options corresponds to scenario 2 of the paper
K0 = 273.15 * Kelvin;
options = struct( ...
    'rateCharge',       6.14062 * 10^5 * meter^3/day/2,   ... % Charge rate (m³/day)
    'rateIdle',         0.0 * kilogram/day/2,             ... % Idle rate (kg/day)
    'rateCushion',      6.14062 * 10^5 * meter^3/day/2,   ... % Cushion (H₂) injection rate (m³/day)
    'rateDischarge',    6.14062 * 10^5 * meter^3/day/2,   ... % Discharge rate (m³/day)
    'bhp',              90.0 * barsa,                     ... % Production bottom hole pressure (BHP) in bar
    'tempCharge',       K0 + 80 * Kelvin,                 ... % Charging temperature (K)
    'tempDischarge',    K0 + 80 * Kelvin,                 ... % Discharging temperature (K)
    'tempCushion',      K0 + 80 * Kelvin,                 ... % Cushion temperature (K)
    'timeCushion',      120 * day,                        ... % Duration of cushioning phase (days)
    'timeCharge',       30 * day,                         ... % Duration of charging phase (days)
    'timeIdle',         15 * day,                         ... % Duration of idle phase (days)
    'timeShut',         30 * day,                         ... % Duration of shut phase (days)
    'timeDischarge',    30 * day,                         ... % Duration of discharging phase (days)
    'dtCharge',         1.0 * day,                        ... % Time step during charging (days)
    'dtCushion',        1.0 * day,                        ... % Time step during cushioning (H₂) (days)
    'dtIdle',           1.0 * day,                        ... % Time step during idle (days)
    'dtShut',           1.0 * day,                        ... % Time step during shut phase (days)
    'dtDischarge',      1.0 * day,                        ... % Time step during discharging (days)
    'numCycles',        20,                               ... % Total number of cycles (charging and discharging)
    'numCyclesCushions',6,                                ... % Number of cycles for cushion gas (H₂)
    'chargeOnly',       0,                                ... % Flag to simulate only the charging phase (0: No, 1: Yes)
    'cushionOnly',      0,                                ... % Flag to simulate only the cushioning phase (0: No, 1: Yes)
    'dischargeOnly',    0,                                ... % Flag to simulate only the discharging phase (0: No, 1: Yes)
    'useGroupCtrl',     false,                            ... % Flag to use group control (true/false)
    'initPres',         81.6 * barsa,                     ... % Initial pressure (bar)
    'initSat',          [0.1,0.1, 0.7999],                           ... % Initial saturation
    'initTemp',         273.15 + 80,                      ... % Initial temperature (°C)
    'use_bc',           true,                             ... % Flag to use boundary conditions (true/false)
    'use_cushion',      true,                             ... % Flag to use cushion phase (true/false)
    'nbact0',           10^7, ....
    'use_bhp',          true,                              ... % Flag to use bottom hole pressure (true/false)
    'initComp',  [0.0045  0.8723 0.0904  0.0250 0.0072 3.9800e-04 1.9900e-04 0]...% [0.8,0.0,0.006,0.018,0.176] ...
    );

%% Process optinal input arguments
[options, fullSetup, ~] = processTestCaseInput(mfilename, ...
    options, description, varargin{:});
options = checkOptions(options);
if ~fullSetup, return; end
%% Optional input arguments
options = merge_options(options, varargin{:});
if nargout <= 2, return;
end
% just to be sure on the module dependencies
require ad-core ad-props ad-blackoil
%% perm, poro, and grid as in the benchmark
G = processGRDECL(deck.GRID, 'checkgrid', false);
G = computeGeometry(G);
deck = convertDeckUnits(deck);
PERMX = deck.GRID.PERMX;
PERMY = deck.GRID.PERMY;
PERMZ = deck.GRID.PERMZ;
perm =  [PERMX PERMY PERMZ];
poro =  deck.GRID.PORO; clear p
rock = makeRock(G, perm, poro);

%% Set up schedule and initiate model from deck
%schedule = setUpSchedule(G, rock, options);
[state0, model, schedule, nls] = initEclipseProblemAD(deck,'getSchedule',true,'getInitialState', false);
deckfluid = model.fluid;
krpts.o = [0.2000 0.7800 0.7800 0];
krpts.g = [0 0.8000 1 0];
% fluid = struct(...
%     'krO', fluid.krW, ...              % Water relative permeability
%     'krG', fluid.krOW, ...                   % Gas relative permeability
%     'krPts', krpts, ...                                    % Placeholder for additional permeability points if needed
%     'cW', fluid.cW,                                         % Compressibility of water (example)
%     'muWr', fluid.muW,                                       % Water viscosity (example)
%     'bW', fluid.bW,                                          % Water formation volume factor (example)
%     'muW', fluid.,                   % Water viscosity dependence on pressure
%     'rhoOS', 0,                                               % Oil density, if relevant (set to zero for gas)
%     'rhoWS', 1.0203e+03,                                      % Water density (kg/m³)
%     'rhoGS', 0.8,                                             % Gas density (example value, can vary with pressure/temperature)
%     'pvMultR', @(p) pvMult(p, cR, pRef)                       % Multiplier function for gas, e.g., pv multiplication with reference
% );
deck = model.inputdata;
% for i=1:length(schedule.control)
%     schedule.control(i).W.compi=[0, 1];
% end
rock.regions.saturation=deck.REGIONS.SATNUM;
gravity reset on;
%% We reset a compositional model and "oil" is the water phase
% compFluid =TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide',...
%     'Nitrogen','Methane'}, {'H2O','H2','CO2','N2','CH4'});

% compFluid = TableCompositionalMixture(...
%     {'Water','Hydrogen','CarbonDioxide','Nitrogen','Methane','Ethane','Propane','Butane'}, ...
%     {'H2O','H2','CO2','N2','CH4','C2H6','C3H8','C4H10'});
compFluid = TableCompositionalMixture(...
    {'Water','Methane','Nitrogen','CarbonDioxide','Ethane','Propane','Butane','Hydrogen'}, ...
    {'H2O','C1','N2','CO2','C2','C3','NC4','H2'});
 
[rhow,rhog]=deal(1000* kilogram/meter^3,8.1688* kilogram/meter^3); %density kilogram/meter^3;
[viscow,viscog]=deal(1.0*centi*poise,0.0094234*centi*poise);%viscosity
[cfw,cfg]=deal(0,8.1533e-3/barsa); %compressibility
 
[srw,src]=deal(0.0,0.0);
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
[Pmaxz,Pref1,Pminz,Pe]=deal(95*barsa,114*barsa,120*barsa,0.1*barsa); %pressions

%initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
% fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
%                          'rho',[rhow,rhog],'pRef',options.initPres,...
%                          'c',[cfw,cfg],'n',[2,2],'smin',[srw,src]);

% 
% arg = {model.G, rock, fluid, compFluid,...
%     'water', false, 'oil', true, 'gas', true,... % water-oil system
% 	'bacteriamodel', true,'diffusioneffect',false,'liquidPhase', 'O',...
%     'vaporPhase', 'G'}; % water=liquid, gas=vapor
% model = BioChemsitryModel(arg{:});
%% Set up initial state
% state0 = setUpInitialState(model, schedule.control(1).W, options);
T0 = options.initTemp;
s0= options.initSat; %initial saturations  Sw=1
z0 = options.initComp; %initial composition: H2O,H2,CO2,N2,CH4.

%===Bacteria model===========================
% if model.bacteriamodel
%     nbact0=options.nbact0.*0;
%     state0 = initCompositionalStateBacteria(model,options.initPres.*ones(G.cells.num,1) ,T0,s0,z0,nbact0);
% else
    state0 = initCompositionalState(model, options.initPres.*ones(G.cells.num,1) , T0, s0, z0);
% end
%% Plotting
plotOptions = {'View'              , [0,0]         , ...
    'PlotBoxAspectRatio', [1,1,0.25]    , ...
    'Projection'        , 'orthographic', ...
    'Size'              , [800, 300]    };
% deck.RUNSPEC.TITLE ='H2_illustration_storage';
% deck_new = model2Deck(model, schedule, 'deck', deck);
end

function W = setUpWells(G, rock, options)

    %% We also reset the well in the highest point: we reset Well coordinates (wc) and radii (r)
    % ps this is slightly different from the benchmark
    wc = [5891; 9612; 13333; 17054; 20775; 24496; 28217; 31938; 35659; 39380; 43101];
    r = [0.9479; 0.7354; 0.7338; 2.2337; 2.2337; 2.6778; 2.6722; 5.1490; 5.1491; 0.7834; 0.7834];
    
    %% Add production wells with specified parameters
    W = addWell([], G, rock, wc, 'Name', 'Prod', ...
                'Radius', r, 'type', 'rate', ...
                'val', options.rateCharge, 'compi', [0, 1]);
   W(1).components = [0.0 0.95 0.05 0.0 0.0]; 
    
    %% Set groups if group control is used
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'});
    end

end

function schedule = setUpSchedule(G0, rock, options)
    %% Set up initial wells based on provided grid, rock, fluid, and options
    W = setUpWells(G0, rock, options);     
    
    %% If cushion stage is being used
    if options.use_cushion
        if options.use_bhp
            % Set well properties for BHP control during cushion phase
            W(1).type     = 'bhp';
            W(1).name     = 'cushion';
            W(1).val      = options.bhp;
            W(1).T        = options.tempCushion;
            W(1).sign     = 1;

            %% Create time steps for cushion phase and set BHP value
            dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
            bhpCushion = options.bhp;

            %% Schedule for the initial 9 cushion steps
            for i = 1:9    
                dtCushion = dtCushions(i);
                W(1).val = bhpCushion;  % Fixed BHP value
                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
            end

            %% Intermediate cushion steps
            dtCushion = dtCushions(10:end-10);
            W(1).val = bhpCushion; 
            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

            %% Final cushion steps
            for i = 1:9    
                dtCushion = dtCushions(end-9+i);
                W(1).val = bhpCushion;  % Fixed BHP value
                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
            end
        else
            %% Set well properties for rate control during cushion phase
            W(1).type     = 'rate';
            W(1).name     = 'cushion';
            W(1).val      = options.rateCushion;
            W(1).T        = options.tempCushion;
            W(1).sign     = 1;

            %% Create time steps and rates for cushion phase
            dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
            rateCushion = options.rateCushion .* dtCushions ./ max(dtCushions);

            %% Schedule for the initial 9 cushion steps with rate control
            for i = 1:9
                dtCushion = dtCushions(i);
                W(1).val = rateCushion(10);
                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
            end

            %% Intermediate cushion steps with constant rate
            dtCushion = dtCushions(10:end-10);
            W(1).val = rateCushion(10); 
            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

            %% Final cushion steps with decreasing rates
            for i = 1:9
                dtCushion = dtCushions(end-9+i);
                W(1).val = rateCushion(end-9+i);  % Vary rate for last steps
                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
            end
        end
    end

    %% Set up wells for charge phase
    W = setUpWells(G0, rock, options);
    W(1).type = 'rate';
    W(1).name = 'charge';
    W(1).val  = options.rateCharge;
    W(1).T    = options.tempCharge;
    W(1).sign = 1;

    %% Create time steps and rates for charge phase
    dtCharges = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
    rateCharge = options.rateCharge .* dtCharges ./ max(dtCharges);

    %% Schedule for the charge phase (initial 9 steps)
    for i = 1:9
        dtCharge = dtCharges(i);
        W(1).val = rateCharge(10);
        scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
    end

    %% Intermediate charge steps with fixed rate
    dtCharge = dtCharges(10:end-10);
    W(1).val = rateCharge(10);
    scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);

    %% Final charge steps with constant rate
    for i = 1:9
        dtCharge = dtCharges(end-9+i);
        W(1).val = rateCharge(10);
        scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
    end

    %% Handle group control if enabled
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end

    %% Set well properties for idle (shut) phase
    W(1).type = 'rate';
    W(1).val  = options.rateIdle;
    W(1).name = 'shut';
    W(1).T    = options.tempCushion;
    W(1).sign = 1;

    %% Create time steps for idle phase
    dtIdle = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
    scheduleIdle = simpleSchedule(dtIdle, 'W', W);
    if options.useGroupCtrl
        groups = [];
        scheduleIdle.groups = groups;
    end

    %% Set up the shut phase
    dtShut = rampupTimestepsEnds(options.timeShut, options.dtShut);
    scheduleShut = simpleSchedule(dtShut, 'W', W);

    %% Set well properties for discharge phase
    W(1).type = 'rate';
    W(1).name = 'discharge';
    W(1).val  = -options.rateDischarge;
    W(1).sign = -1;
    W(1).lims.bhp = 20*barsa;
    W(1).cstatus(2:end) = 0;

    %% Create time steps for discharge phase
    dtDischarge = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);

    %% Handle group control for discharge
    if options.useGroupCtrl
        groups = [];
        scheduleCharge.groups = groups;
    end

    %% Combine schedules based on user options (charge, discharge, or cushion-only)
    if options.chargeOnly    
        schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    elseif (options.cushionOnly && false)
        schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
    else
        if options.use_cushion
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);

            %% Set different cushion values and combine
            scheduleCushion1 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 94*barsa();
            end
            scheduleCushion2 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 98*barsa();
            end
            scheduleCushion3 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
            for i = 1:length(scheduleCushions)
                scheduleCushions{i}.control.W.val = 102*barsa();
            end
            scheduleCushion4 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);

            %% Final combination of schedules
            if options.cushionOnly
                schedule = combineSchedules(scheduleCushion1, scheduleCushion2, scheduleCushion3, scheduleCushion4, 'makeConsistent', false);
            else
                schedule = combineSchedules(scheduleCushion1, scheduleCushion2, scheduleCushion3, scheduleCushion4, scheduleCushion4, scheduleCushion4, schedule{:}, 'makeConsistent', false);
            end    
        else
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
        end
    end

    %% Apply boundary conditions if enabled
    if options.use_bc
        bc = setUpBc(G0, options);
        for i = 1:numel(schedule.control)
            schedule.control(i).bc = bc;
        end
    end
end

function bc = setUpBc(G, options)
    %% Get boundary faces of the grid
    f = boundaryFaces(G);
    
    %% Select lateral faces
    f1 = (G.faces.normals(f, 3) > 2.48);
    faces = f(~f1);

    %% Calculate the displacement from a reference point (1140 units)
    dx = bsxfun(@minus, G.faces.centroids(faces, :), 1140);
    
    %% Calculate pressure drop due to gravity using omega and displacement
    rhoOS = 9.98e+02;
    dp = rhoOS .* (dx * reshape(gravity, [], 1));
    
    %% Define the pressure at the boundary (init pressure plus the gravity-induced pressure drop)
    pressure = options.initPres + dp;

    %% Set the boundary condition with the computed pressure and saturation
    bc = addBC([], faces, 'pressure', pressure, 'sat', options.initSat);   
    bc.components = repmat(options.initComp, numel(bc.face), 1);
end



% function schedule = setUpSchedule(G0, rock, options)
%         
%     W = setUpWells(G0, rock, options);     
%     if options.use_cushion
%        if options.use_bhp
%           W(1).type     = 'bhp';
%           W(1).name     = 'cushion';
%           W(1).val      = options.bhp;
%           W(1).T        = options.tempCushion;
%           W(1).sign     = 1;
%     
%           dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%           bhpCushion = [options.bhp];
%           for i = 1:9    
%               dtCushion      = dtCushions(i);
%               W(1).val      = bhpCushion(1);
%               scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%           end
% 
%           dtCushion       = dtCushions(10:end-10);
%           W(1).val      = bhpCushion(end);
%           scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%           for i = 1:9    
%               dtCushion      = dtCushions(end-9+i);
%               W(1).val      = bhpCushion(end);
%               scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%           end
%        else
%            
%            W(1).type     = 'rate';
%            W(1).name     = 'cushion';
%            W(1).val      = options.rateCushion;
%            W(1).T        = options.tempCushion;
%            W(1).sign     = 1;
%     
%            dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%            rateCushion = options.rateCushion.*dtCushions./max(dtCushions);
%            for i = 1:9               
%                dtCushion      = dtCushions(i);
%                W(1).val      = rateCushion(10);
%                scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%            end
% 
%            dtCushion       = dtCushions(10:end-10);
%            W(1).val      = rateCushion(10);
%            scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%        
%            for i = 1:9               
%                dtCushion      = dtCushions(end-9+i);
%                W(1).val      = rateCushion(end-9+i);
%                scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%            end
%        end
%     end
%    
%     W = setUpWells(G0, rock, options);    
%     W(1).type     = 'rate';
%     W(1).name     = 'charge';    
%     W(1).val      = options.rateCharge;
%     W(1).T        = options.tempCharge;
%     W(1).sign     = 1;
% %     W.lims.rate = options.rateCharge;
%     dtCharges       = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
%     rateCharge = options.rateCharge.*dtCharges./max(dtCharges);
%     for i = 1:9    
%         dtCharge      = dtCharges(i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
%     end
% 
%     dtCharge       = dtCharges(10:end-10);
%     W(1).val      = rateCharge(10);
%     scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);
% 
%        
%     for i = 1:9    
%         dtCharge      = dtCharges(end-9+i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
%     end
%     
%     
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
% 
%     W(1).type     = 'rate';
%     W(1).val      = options.rateIdle;
%     W(1).name     = 'shut';        
%     W(1).T        = options.tempCushion;
%     W(1).sign     = 1;
%     
%     dtIdle       = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
%     scheduleIdle = simpleSchedule(dtIdle, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleIdle.groups = groups;
%     end
% 
%     dtShut       = rampupTimestepsEnds(options.timeShut, options.dtShut);
%     scheduleShut = simpleSchedule(dtShut, 'W', W);
% 
%     W(1).type     = 'rate';
%     W(1).name     = 'discharge';    
%     W(1).val      = -options.rateDischarge;
%     W(1).sign     = -1;
%     W(1).lims.bhp = 20*barsa;
%     W(1).cstatus(2:end) = 0;
%     dtDischarge       = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
%     scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
%     
%     if options.chargeOnly    
%         schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
%     elseif options.dischargeOnly
%         schedule = scheduleDischarge;
%     elseif (options.cushionOnly&&false)
%         schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
%     else 
%         if options.use_cushion
%            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%            schedule = repmat({schedule}, 1, options.numCycles);
%            scheduleCushion1 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
%            for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =94*barsa();
%            end
%            scheduleCushion2 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false); 
%             for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =98*barsa();
%            end
%            scheduleCushion3 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);
%            for i =1:length(scheduleCushions)
%               scheduleCushions{i}.control.W.val =102*barsa();
%            end
%            scheduleCushion4 = combineSchedules(scheduleCushions{:}, scheduleShut, 'makeConsistent', false);  
%            if options.cushionOnly           
%                schedule = combineSchedules(scheduleCushion1,scheduleCushion2,scheduleCushion3,scheduleCushion4, 'makeConsistent', false);
%            else              
%                schedule = combineSchedules(scheduleCushion1,scheduleCushion2,scheduleCushion3,scheduleCushion4, scheduleCushion4, scheduleCushion4, schedule{:}, 'makeConsistent', false);
%            end    
%         else
%             schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%             schedule = repmat({schedule}, 1, options.numCycles);
%             schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
%         end
%     end
% 
%         
%     if options.use_bc    
%         bc = setUpBc(G0,rock,fluid,options);        
%         for i = 1:numel(schedule.control)        
%             schedule.control(i).bc = bc;
%         end
%     end
% 
% end


%-------------------------------------------------------------------------%
function state0 = setUpInitialState(model, W, options)
    %% Initialize the reservoir solution based on the model grid and initial pressure
    state0 = initResSol(model.G, options.initPres, options.initSat);    
    state0.rs = zeros(size(state0.pressure));  % Set initial solution gas-to-oil ratio to zero
    %% Initialize well solutions
    wellSol = initWellSolAD(W, model, state0);
    wellSol.bhp = options.initPres;  % Set initial bottom hole pressure for all wells
    state0.wellSol = wellSol;  % Store well solutions in the state structure
end


%-------------------------------------------------------------------------%
function options = checkOptions(options)
    
    assert(~(options.chargeOnly && options.dischargeOnly), ...
        'Cannot simulate only charge and only discharge at the same time');
    
end

function dT = rampupTimestepsEnds(time, dt, n)
% Create timesteps that ramp up geometrically
%
% SYNOPSIS:
%   dT = rampupTimesteps(1*year, 30*day)
%   dT = rampupTimesteps(1*year, 30*day, 5)
%
% DESCRIPTION:
%   This function generates a timestep sequence for a given total time
%   interval that increases geometrically until it reaches some target
%   timestep. The rest of the interval is then divided into a number of
%   target timesteps.
%
% REQUIRED PARAMETERS:
%   time   - The total simulation time so that sum(dt) = time
%
%   dt     - Target timestep after initial ramp-up
%
%   n      - (OPTIONAL) Number of rampup steps. Defaults to 8.
%
% RETURNS:
%   dt     - Array of timesteps.
%
% NOTE:
%   The final timestep may be shorter than dt in order to exactly reach T.
%

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if nargin < 3
        n = 1;
    end
    if time == 0
        dT = [];
        return
    end
    % Initial geometric series
    dt_init = (dt./2.^[n n:-1:1])';
    cs_time = cumsum(dt_init);
    if any(cs_time > time)
        dt_init = dt_init(cs_time < time);
    end
    
    % Remaining time that must be discretized
    dt_left = time - sum(2.*dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(2.*dt_init) - sum(dt_rem);
    % Less than to account for rounding errors leading to a very small
    % negative time-step.
    if dt_final <= 0
        dt_final = [];
    end
       
    if dt_final >= dt_init(1)
        dt_final = dt_init(1);
    end
    % Combined timesteps
    dT = [dt_init; dt_rem;sort(dt_init,'descend'); dt_final];
end