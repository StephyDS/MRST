%% MRST Simulation for Hydrogen Storage with Bacterial Growth
%
% Description:
% This MRST example simulates underground hydrogen storage with and without
% microbial activity in a 3D porous medium using compositional fluid
% properties. The case is based on:
%
% Khoshnevis, N., Hogeweg, S., Goncalves Machado, C., and Hagemann, B.
% "Numerical Modeling of Bio-Reactive Transport During Underground
% Hydrogen Storage – a Benchmark Study", EAGE, vol. 1, pp. 1-5, 2023.
% https://doi.org/10.3997/2214-4609.202321087
%
% Features:
% - Two-phase system (liquid 'O' and gas 'G')
% - Four components (H2O, H2, CO2, CH4)
% - With/without microbial activity of archaea
% - Cyclic injection/production schedule
% - Based on EAGE 2023 Benchmark case
%
% Author: Stéphanie Delage Santacreu
% Date: 23/09/2025
% Organization: Université de Pau et des Pays de l'Adour, E2S UPPA,
%               CNRS, LFCR, UMR5150, Pau, France
% -------------------------------------------------------------------------

%% Initialization
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
mrstVerbose = true;
gravity reset on

%% ============ Grid and Rock Properties =====================
[nx, ny, nz] = deal(15, 15, 8);
[Lx, Ly, Lz] = deal(1525, 1525, 50);

G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 1000; % Reservoir depth
G = computeGeometry(G);

perm = [100, 100, 10] .* milli*darcy;
rock = makeRock(G, perm, 0.2);

%% Fluid Properties
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'C1'});

[rhow, rhog]   = deal(999.7 * kilogram/meter^3, 1.2243 * kilogram/meter^3);
[viscow, viscog] = deal(1.3059 * centi*poise, 0.01763 * centi*poise);
[cfw, cfg]     = deal(5.0e-5/barsa, 1.0/barsa);

[srw, src] = deal(0.2, 0.05);
P0 = 100 * barsa;

fluid = initSimpleADIFluid('phases', 'OG', ...
    'mu', [viscow, viscog], ...
    'rho', [rhow, rhog], ...
    'pRef', P0, ...
    'c', [cfw, cfg], ...
    'n', [2, 2], ...
    'smin', [srw, src]);

Pe = 0.1 * barsa;
fluid.pcOG = @(sg) Pe * max((1 - sg - srw) ./ (1 - srw), 1e-5).^(-1/2);

%% Wells
W = [];
n1 = floor(0.5*nx) + 1;
n2 = floor(0.5*ny) + 1;

% Injector (CO2-rich)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'comp_i', [0, 1], 'Radius', 0.5, 'name', 'Injector', ...
    'type', 'rate', 'Val', 1e6*meter^3/day, 'sign', 1);
W(end).components = [0.0, 0.6, 0.4, 0.0];

% Rest (shut-in)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'compi', [0, 1], 'Radius', 0.5, 'name', 'Rest', ...
    'type', 'rate', 'Val', 0.0, 'sign', 1);
W(end).components = [0.0, 0.95, 0.05, 0.0];

% Injector (H2-rich)
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'comp_i', [0, 1], 'Radius', 0.5, 'name', 'InjectorH2', ...
    'type', 'rate', 'Val', 1e6*meter^3/day, 'sign', 1);
W(end).components = [0.0, 0.95, 0.05, 0.0];

% Producer
W = verticalWell(W, G, rock, n1, n2, 1, ...
    'compi', [0, 1], 'Radius', 0.5, 'name', 'Producer', ...
    'type', 'rate', 'Val', -1e6*meter^3/day, 'sign', -1);
W(end).components = [0.0, 0.95, 0.05, 0.0];

%% Time Schedule
ncycles    = 3; 
deltaT     = 1*day; %5*day;
nbj_buildUp = 60*day; nbj_rest = 20*day;
nbj_inject  = 30*day; nbj_idle = 20*day;
nbj_prod    = 30*day; nbj_idle1 = 20*day;

[schedule, TotalTime, nbuildUp, nrest, ninject, nidle, nprod, nidle1] = ...
    createCyclicScenario(deltaT, ncycles, nbj_buildUp, nbj_rest, ...
    nbj_inject, nbj_idle, nbj_prod, nbj_idle1, W);
% smaller step for initialization
schedule.step.val(1) = 1*day;

%% EOS and Model
eosname = 'sw';
compEOS = EquationOfStateModel(G, compFluid, eosname);

backend = DiagonalAutoDiffBackend('modifyOperators', true);
%mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);

%% Initial State
T0 = 40 + 273.15;
s0 = [0.2, 0.8];
z0 = [0.7, 0.0, 0.02, 0.28];

%% --- Simulation 0: Without bacteria no molecular diffusion---
arg = {G, rock, fluid, compFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', false, 'moleculardiffusion',false, 'dispersion',false,...
    'liquidPhase', 'O', 'vaporPhase', 'G'};

model_nobact = BiochemistryModel(arg{:});
model_nobact.outputFluxes = false;
model_nobact.EOSModel = compEOS;

state0_nobact = initCompositionalState(model_nobact, P0, T0, s0, z0);

lsolve = selectLinearSolverAD(model_nobact);
nls = NonLinearSolver(); nls.LinearSolver = lsolve;


problem_nobact = packSimulationProblem(state0_nobact, model_nobact, schedule, ...
    'Benchmark_NoBacteria_6cycles_15_dt1', 'NonLinearSolver', nls);
simulatePackedProblem(problem_nobact)%;, 'restartStep',1);
[ws_nobact, states_nobact] = getPackedSimulatorOutput(problem_nobact);
results_nobact = postProcessResults(states_nobact, ws_nobact, model_nobact, 'nobact');


%% --- Simulation 1: Without bacteria molecular diffusion---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', false, 'moleculardiffusion',true,'dispersion',true,...
%     'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_nobact_moldiff = BiochemistryModel(arg{:});
% model_nobact_moldiff.outputFluxes = false;
% model_nobact_moldiff.EOSModel = compEOS;
% 
% state0_nobact_moldiff = initCompositionalState(model_nobact_moldiff, P0, T0, s0, z0);
% 
% lsolve = selectLinearSolverAD(model_nobact_moldiff);
% nls = NonLinearSolver(); nls.LinearSolver = lsolve;
% model_nobact_moldiff.nonlinearTolerance=1.e-5;
% 
% problem_nobact_moldiff = packSimulationProblem(state0_nobact_moldiff, model_nobact_moldiff, schedule, ...
%     'Benchmark_NoBacteria_moldiff_6cycles_21_dt1', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_nobact_moldiff);%, 'restartStep',1);
% [ws_nobact_moldiff, states_nobact_moldiff] = getPackedSimulatorOutput(problem_nobact_moldiff);
% results_nobact_moldiff = postProcessResults(states_nobact_moldiff, ws_nobact_moldiff, model_nobact_moldiff, 'nobact');
% 

%% --- Simulation 2: Without bacteria salt water---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', false, 'moleculardiffusion',false,'dispersion',false,...
%     'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_nobact_salt = BiochemistryModel(arg{:});
% model_nobact_salt.outputFluxes = false;
% model_nobact_salt.EOSModel = compEOS;
% model_nobact_salt.EOSModel.msalt=3;
% 
% state0_nobact_salt = initCompositionalState(model_nobact_salt, P0, T0, s0, z0);
% 
% lsolve = selectLinearSolverAD(model_nobact_salt);
% nls = NonLinearSolver(); nls.LinearSolver = lsolve;
% 
% problem_nobact_salt = packSimulationProblem(state0_nobact_salt, model_nobact_salt, schedule, ...
%     '3', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_nobact_salt, 'restartStep',1);
% [ws_nobact_salt, states_nobact_salt] = getPackedSimulatorOutput(problem_nobact_salt);
% results_nobact_salt = postProcessResults(states_nobact_salt, ws_nobact_salt, model_nobact_salt, 'nobact');


% %% --- Simulation 2: Without bacteria, dispersion---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', false, 'moleculardiffusion',false,'dispersion',true,...
%     'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_nobact_disp = BiochemistryModel(arg{:});
% model_nobact_disp.outputFluxes = false;
% model_nobact_disp.EOSModel = compEOS;
% 
% state0_nobact_disp = initCompositionalState(model_nobact_disp, P0, T0, s0, z0);
% 
% lsolve = selectLinearSolverAD(model_nobact_disp);
% nls = NonLinearSolver(); nls.LinearSolver = lsolve;
% 
% problem_nobact_disp = packSimulationProblem(state0_nobact_disp, model_nobact_disp, schedule, ...
%     'Benchmark_NoBacteria_disp_6cycles_31_dt1', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_nobact_disp);%, 'restartStep',1);
% [ws_nobact_disp, states_nobact_disp] = getPackedSimulatorOutput(problem_nobact_disp);
% 
% results_nobact_disp = postProcessResults(states_nobact_disp, ws_nobact_disp, model_nobact_disp, 'nobact');
% 

%% --- Simulation 3: Without bacteria, dispersion and moldiff---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', false, 'moleculardiffusion',true,'dispersion',true,...
%     'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_nobact_diffdisp = BiochemistryModel(arg{:});
% model_nobact_diffdisp.outputFluxes = false;
% model_nobact_diffdisp.EOSModel = compEOS;
% 
% state0_nobact_diffdisp = initCompositionalState(model_nobact_diffdisp, P0, T0, s0, z0);
% 
% lsolve = selectLinearSolverAD(model_nobact_diffdisp);
% nls = NonLinearSolver(); nls.LinearSolver = lsolve;
% %nls = NonLinearSolver('useRelaxation', true,'Verbose', true);nls.LinearSolver = lsolve;
% 
% problem_nobact_diffdisp = packSimulationProblem(state0_nobact_diffdisp, model_nobact_diffdisp, schedule, ...
%     'Benchmark_NoBacteria_diffdisp_6cycles_21b', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_nobact_diffdisp);%, 'restartStep',1);
% [ws_nobact_diffdisp, states_nobact_diffdisp] = getPackedSimulatorOutput(problem_nobact_diffdisp);
% 
% results_nobact_diffdisp = postProcessResults(states_nobact_diffdisp, ws_nobact_diffdisp, model_nobact_diffdisp, 'nobact');


%% --- Efficiency and Comparison ---
eff_nobact = calculateH2Efficiency(results_nobact.H2_well, ...
    nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
fprintf('H2 Production Efficiency (No bacteria): %.4f%%\n', eff_nobact(end));
% 
% eff_nobact_moldiff = calculateH2Efficiency(results_nobact_moldiff.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (No bacteria, With molecular diffusion): %.4f%%\n', eff_nobact_moldiff(end));
% 
% eff_nobact_salt = calculateH2Efficiency(results_nobact_salt.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (No bacteria, With salt water): %.4f%%\n', eff_nobact_salt(end));

% eff_nobact_disp = calculateH2Efficiency(results_nobact_disp.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (No bacteria, With dispersion): %.4f%%\n', eff_nobact_disp'end));


% Component changes
% H2_loss_moldiff = (abs(results_nobact.totMassH2 - results_nobact_moldiff.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_moldiff = (abs(results_nobact.totMassCO2 - results_nobact_moldiff.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_moldiff  = (abs(results_nobact.totMassC1 - results_nobact_moldiff.totMassC1) ./ results_nobact.totMassC1) * 100;
% % 
% H2_loss_salt = ((-results_nobact.totMassH2 + results_nobact_salt.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_salt = ((-results_nobact.totMassCO2 + results_nobact_salt.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_salt  = ((-results_nobact.totMassC1 + results_nobact_salt.totMassC1) ./ results_nobact.totMassC1) * 100;

% H2_loss_disp = (abs(results_nobact.totMassH2 - results_nobact_disp.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_disp = (abs(results_nobact.totMassCO2 - results_nobact_disp.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_disp  = (abs(results_nobact.totMassC1 - results_nobact_disp.totMassC1) ./ results_nobact.totMassC1) * 100;
% 
% fprintf('Total H2 loss due to molecular diffusion: %.4f%%\n', H2_loss_moldiff(end));
% fprintf('Total CO2 loss due to molecular diffusion: %.4f%%\n', CO2_loss_moldiff(end));
% fprintf('Total C1 gain due to molecular diffusion:  %.4f%%\n', C1_gain_moldiff(end));

% fprintf('Total H2 loss due to salt water: %.4f%%\n', H2_loss_salt(end));
% fprintf('Total CO2 loss due to salt water: %.4f%%\n', CO2_loss_salt(end));
% fprintf('Total C1 gain due to salt water:  %.4f%%\n', C1_gain_salt(end));
% 
% fprintf('Total H2 loss due to dispersion: %.4f%%\n', H2_loss_disp(end));
% fprintf('Total CO2 loss due to dispersion: %.4f%%\n', CO2_loss_disp(end));
% fprintf('Total C1 gain due to dispersion:  %.4f%%\n', C1_gain_disp(end));



%% --- Simulation 5: With bacteria pure water no diffusion, no dispersion---
arg = {G, rock, fluid, compFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',false,...
    'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};

model_bact = BiochemistryModel(arg{:});

model_bact.outputFluxes = false;
model_bact.EOSModel = compEOS;

nbact0 = 1; model_bact.nbactMax = 1e8;
state0_bact = initCompositionalStateBacteria(model_bact, P0, T0, s0, z0, nbact0, compEOS);

lsolve = selectLinearSolverAD(model_bact);
nls.LinearSolver = lsolve;

problem_bact = packSimulationProblem(state0_bact, model_bact, schedule, ...
    'Benchmark_Bacteria_6cycles_15_dt1', 'NonLinearSolver', nls);
simulatePackedProblem(problem_bact);%, 'restartStep',1);
[ws_bact, states_bact] = getPackedSimulatorOutput(problem_bact);
results_bact = postProcessResults(states_bact, ws_bact, model_bact, 'bact');


%% --- Simulation 6: With bacteria pure water diffusive bact, no diffusion, no dispersion---
arg = {G, rock, fluid, compFluid, true, backend, ...
    'water', false, 'oil', true, 'gas', true, ...
    'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',false,...
    'bactdiffusion',true,'chemotaxisEffect',true,'liquidPhase', 'O', 'vaporPhase', 'G'};

model_bact_diffbact = BiochemistryModel(arg{:});

model_bact_diffbact.outputFluxes = false;
model_bact_diffbact.EOSModel = compEOS;

nbact0 = 1; model_bact_diffbact.nbactMax = 1e8;
state0_bact_diffbact = initCompositionalStateBacteria(model_bact_diffbact, P0, T0, s0, z0, nbact0, compEOS);

lsolve = selectLinearSolverAD(model_bact_diffbact);
nls.LinearSolver = lsolve;

problem_bact_diffbact = packSimulationProblem(state0_bact_diffbact, model_bact_diffbact, schedule, ...
    'Benchmark_Bacteria_chemo_6cycles_15_dt1', 'NonLinearSolver', nls);
simulatePackedProblem(problem_bact_diffbact, 'restartStep',1);
[ws_bact_diffbact, states_bact_diffbact] = getPackedSimulatorOutput(problem_bact_diffbact);
results_bact_diffbact = postProcessResults(states_bact_diffbact, ws_bact_diffbact, model_bact_diffbact, 'bact');

%% --- Simulation 6: With bacteria pure water WITH diffusion, no dispersion---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', true, 'moleculardiffusion',true,'dispersion',true,...
%     'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_bact_moldiff = BiochemistryModel(arg{:});
% 
% model_bact_moldiff.outputFluxes = false;
% model_bact_moldiff.EOSModel = compEOS;
% 
% nbact0 = 1; model_bact_moldiff.nbactMax = 1e8;
% state0_bact_moldiff = initCompositionalStateBacteria(model_bact_moldiff, P0, T0, s0, z0, nbact0, compEOS);
% 
% lsolve = selectLinearSolverAD(model_bact_moldiff);
% nls.LinearSolver = lsolve;
% model_bact.nonlinearTolerance=1.e-5;
% 
% problem_bact_moldiff = packSimulationProblem(state0_bact_moldiff, model_bact_moldiff, schedule, ...
%     'Benchmark_Bacteria_moldiff_6cycles_21_dt1', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_bact_moldiff, 'restartStep',1);
% [ws_bact_moldiff, states_bact_moldiff] = getPackedSimulatorOutput(problem_bact_moldiff);
% results_bact_moldiff = postProcessResults(states_bact_moldiff, ws_bact_moldiff, model_bact_moldiff, 'bact');
% 

% %% --- Simulation 7: With bacteria salt water no diffusion, no dispersion---
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',false,...
%     'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_bact_salt = BiochemistryModel(arg{:});
% 
% model_bact_salt.outputFluxes = false;
% model_bact_salt.EOSModel = compEOS;
% model_nobact_salt.EOSModel.msalt=3;
% 
% nbact0 = 1; model_bact_salt.nbactMax = 1e8;
% state0_bact_salt = initCompositionalStateBacteria(model_bact_salt, P0, T0, s0, z0, nbact0, compEOS);
% 
% lsolve = selectLinearSolverAD(model_bact_salt);
% nls.LinearSolver = lsolve;
% 
% problem_bact_salt = packSimulationProblem(state0_bact_salt, model_bact_salt, schedule, ...
%     'Benchmark_Bacteria_salt3_6cycles_31_dt1', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_bact_salt);%, 'restartStep',1);
% [ws_bact_salt, states_bact_salt] = getPackedSimulatorOutput(problem_bact_salt);
% results_bact_salt = postProcessResults(states_bact_salt, ws_bact_salt, model_bact_salt, 'bact');
% 

%% --- Simulation 8: With bacteria pure water no diffusion, no dispersion-n0=1.e9--
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',false,...
%     'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_bact9 = BiochemistryModel(arg{:});
% 
% model_bact9.outputFluxes = false;
% model_bact9.EOSModel = compEOS;
% 
% nbact0 = 1; model_bact9.nbactMax = 1e9;
% state0_bact9 = initCompositionalStateBacteria(model_bact9, P0, T0, s0, z0, nbact0, compEOS);
% 
% lsolve = selectLinearSolverAD(model_bact9);
% nls.LinearSolver = lsolve;
% 
% problem_bact9 = packSimulationProblem(state0_bact9, model_bact9, schedule, ...
%     'Benchmark_Bacteria_salt3_6cycles_31_dt1_n09', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_bact9);%, 'restartStep',1);
% [ws_bact9, states_bact9] = getPackedSimulatorOutput(problem_bact9);
% results_bact9 = postProcessResults(states_bact9, ws_bact9, model_bact9, 'bact');

%% --- Simulation 9: With bacteria salt water no diffusion, no dispersion-n0=1.e9--
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',false,...
%     'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_bact9_salt = BiochemistryModel(arg{:});
% 
% model_bact9_salt.outputFluxes = false;
% model_bact9_salt.EOSModel = compEOS;
% model_bact9_salt.EOSModel.msalt=3;
% 
% nbact0 = 1; model_bact9_salt.nbactMax = 1e9;
% state0_bact9_salt = initCompositionalStateBacteria(model_bact9_salt, P0, T0, s0, z0, nbact0, compEOS);
% 
% lsolve = selectLinearSolverAD(model_bact9_salt);
% nls.LinearSolver = lsolve;
% 
% problem_bact9_salt = packSimulationProblem(state0_bact9_salt, model_bact9_salt, schedule, ...
%     'Benchmark_Bacteria_salt3_6cycles_31_dt1_n09', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_bact9_salt);%, 'restartStep',1);
% [ws_bact9_salt, states_bact9_salt] = getPackedSimulatorOutput(problem_bact9_salt);
% results_bact9_salt = postProcessResults(states_bact9_salt, ws_bact9_salt, model_bact9_salt, 'bact');
% 





%% --- Simulation 9: With bacteria pure water no diffusion, no dispersion-n0=1.e8--
% arg = {G, rock, fluid, compFluid, true, backend, ...
%     'water', false, 'oil', true, 'gas', true, ...
%     'bacteriamodel', true, 'moleculardiffusion',false,'dispersion',true,...
%     'bactdiffusion',false,'liquidPhase', 'O', 'vaporPhase', 'G'};
% 
% model_bact_disp = BiochemistryModel(arg{:});
% 
% model_bact_disp.outputFluxes = false;
% model_bact_disp.EOSModel = compEOS;
% 
% nbact0 = 1; model_bact_disp.nbactMax = 1e8;
% state0_bact_disp = initCompositionalStateBacteria(model_bact_disp, P0, T0, s0, z0, nbact0, compEOS);
% 
% lsolve = selectLinearSolverAD(model_bact_disp);
% nls.LinearSolver = lsolve;
% 
% problem_bact_disp = packSimulationProblem(state0_bact_disp, model_bact_disp, schedule, ...
%     'Benchmark_Bacteria_disp_6cycles_31_dt1', 'NonLinearSolver', nls);
% simulatePackedProblem(problem_bact_disp);%, 'restartStep',1);
% [ws_bact_disp, states_bact_disp] = getPackedSimulatorOutput(problem_bact_disp);
% results_bact_disp = postProcessResults(states_bact_disp, ws_bact_disp, model_bact_disp, 'bact');
% 

%% --- Efficiency and Comparison ---
eff_bact = calculateH2Efficiency(results_bact.H2_well, ...
    nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
fprintf('H2 Production Efficiency (with bacteria): %.4f%%\n', eff_bact(end));
eff_bact_diffbact = calculateH2Efficiency(results_bact_diffbact.H2_well, ...
    nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
fprintf('H2 Production Efficiency (with bacteria): %.4f%%\n', eff_bact_diffbact(end));
% eff_bact_moldiff = calculateH2Efficiency(results_bact_moldiff.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (with bacteria+diffusion): %.4f%%\n', eff_bact_moldiff(end));
%
% eff_bact_salt = calculateH2Efficiency(results_bact_salt.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (with bacteria+salt): %.4f%%\n', eff_bact_salt(end));
% 
% eff_bact9 = calculateH2Efficiency(results_bact9.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (with bacteria 1e9): %.4f%%\n', eff_bact9(end));
% 
% eff_bact9_salt = calculateH2Efficiency(results_bact9_salt.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (with bacteria 1e9 and salt): %.4f%%\n', eff_bact9_salt(end));
% eff_bact_disp = calculateH2Efficiency(results_bact_disp.H2_well, ...
%     nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles);
% fprintf('H2 Production Efficiency (with bacteria 1e9): %.4f%%\n', eff_bact_disp(end));

% 
H2_loss_bact = (abs(-results_nobact.totMassH2 + results_bact.totMassH2) ./ results_nobact.totMassH2) * 100;
CO2_loss_bact = (abs(-results_nobact.totMassCO2 + results_bact.totMassCO2) ./ results_nobact.totMassCO2) * 100;
C1_gain_bact  = (abs(-results_nobact.totMassC1 + results_bact.totMassC1) ./ results_nobact.totMassC1) * 100;
H2_loss_bact_diffbact = (abs(-results_nobact.totMassH2 + results_bact_diffbact.totMassH2) ./ results_nobact.totMassH2) * 100;
CO2_loss_bact_diffbact = (abs(-results_nobact.totMassCO2 + results_bact_diffbact.totMassCO2) ./ results_nobact.totMassCO2) * 100;
C1_gain_bact_diffbact  = (abs(-results_nobact.totMassC1 + results_bact_diffbact.totMassC1) ./ results_nobact.totMassC1) * 100;
% C1_gain_bact2  = ((-results_nobact.totMassC1 + results_bact.totMassC1) ./ results_nobact.totMassC1) * 100;
% H2_loss_bact = (abs(results_nobact.totMassH2 - results_bact.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact = (abs(results_nobact.totMassCO2 - results_bact.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact  = (abs(results_nobact.totMassC1 - results_bact.totMassC1) ./ results_nobact.totMassC1) * 100;
% 
% H2_loss_bact_salt = (abs(-results_nobact.totMassH2 + results_bact_salt.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact_salt = (abs(-results_nobact.totMassCO2 + results_bact_salt.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact_salt  = (abs(-results_nobact.totMassC1 + results_bact_salt.totMassC1) ./ results_nobact.totMassC1) * 100;
% % 
% H2_loss_bact9 = (abs(-results_nobact.totMassH2 + results_bact9.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact9 = (abs(-results_nobact.totMassCO2 + results_bact9.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact9  = (abs(-results_nobact.totMassC1 + results_bact9.totMassC1) ./ results_nobact.totMassC1) * 100;
% 
% H2_loss_bact9_salt = (abs(-results_nobact.totMassH2 + results_bact9_salt.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact9_salt = (abs(-results_nobact.totMassCO2 + results_bact9_salt.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact9_salt  = (abs(-results_nobact.totMassC1 + results_bact9_salt.totMassC1) ./ results_nobact.totMassC1) * 100;
% 
% H2_loss_bact_disp = (abs(results_nobact.totMassH2 - results_bact_disp.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact_disp = (abs(results_nobact.totMassCO2 - results_bact_disp.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact_disp  = (abs(results_nobact.totMassC1 - results_bact_disp.totMassC1) ./ results_nobact.totMassC1) * 100;
% % 
% H2_loss_bact_moldiff = (abs(results_nobact.totMassH2 - results_bact_moldiff.totMassH2) ./ results_nobact.totMassH2) * 100;
% CO2_loss_bact_moldiff = (abs(results_nobact.totMassCO2 - results_bact_moldiff.totMassCO2) ./ results_nobact.totMassCO2) * 100;
% C1_gain_bact_moldiff  = (abs(results_nobact.totMassC1 - results_bact_moldiff.totMassC1) ./ results_nobact.totMassC1) * 100;
% 

fprintf('Total H2 loss due to microbes, no diffusion: %.4f%%\n', H2_loss_bact(end));
fprintf('Total CO2 loss due to microbes, no diffusion: %.4f%%\n', CO2_loss_bact(end));
fprintf('Total C1 gain due to microbes, no diffusion:  %.4f%%\n', C1_gain_bact(end));
fprintf('Total H2 loss due to microbes,  diffbact: %.4f%%\n', H2_loss_bact_diffbact(end));
fprintf('Total CO2 loss due to microbes, diffbact: %.4f%%\n', CO2_loss_bact_diffbact(end));
fprintf('Total C1 gain due to microbes, diffbact:  %.4f%%\n', C1_gain_bact_diffbact(end));
% 
% 
% fprintf('Total H2 loss due to microbes + diffusion: %.4f%%\n', H2_loss_bact_moldiff(end));
% fprintf('Total CO2 loss due to microbes+diffusion: %.4f%%\n', CO2_loss_bact_moldiff(end));
% fprintf('Total C1 gain due to microbes+ diffusion:  %.4f%%\n', C1_gain_bact_moldiff(end));

% 
% fprintf('Total H2 loss due to microbes + salt: %.4f%%\n', H2_loss_bact_salt(end));
% fprintf('Total CO2 loss due to microbes+salt: %.4f%%\n', CO2_loss_bact_salt(end));
% fprintf('Total C1 gain due to microbes+ salt:  %.4f%%\n', C1_gain_bact_salt(end));
% % 
% fprintf('Total H2 loss due to microbes 1e9: %.4f%%\n', H2_loss_bact9(end));
% fprintf('Total CO2 loss due to microbes 1e9: %.4f%%\n', CO2_loss_bact9(end));
% fprintf('Total C1 gain due to microbes 1e9:  %.4f%%\n', C1_gain_bact9(end));
% % 
% fprintf('Total H2 loss due to microbes 1e9 and salt: %.4f%%\n', H2_loss_bact9_salt(end));
% fprintf('Total CO2 loss due to microbes 1e9and salt: %.4f%%\n', CO2_loss_bact9_salt(end));
% fprintf('Total C1 gain due to microbes 1e9and salt:  %.4f%%\n', C1_gain_bact9_salt(end));
% 
% fprintf('Total H2 loss due to microbes 1e9: %.4f%%\n', H2_loss_bact_disp(end));
% fprintf('Total CO2 loss due to microbes 1e9: %.4f%%\n', CO2_loss_bact_disp(end));
% fprintf('Total C1 gain due to microbes 1e9:  %.4f%%\n', C1_gain_bact_disp(end));
% 
% 
% f20=figure('Name','Efficiency of H2 production','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:ncycles,eff_nobact,'r','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:ncycles,eff_bact,'b','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:ncycles,eff_bact_moldiff,'k--','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:ncycles,eff_bact9,'m:','MarkerSize',11,'LineWidth',2)
% title(' Efficiency of H2 production','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'cycles'},'FontWeight','bold','Color','k')
% ylabel({'Efficiency (%)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.XMinorTick='on';
% ax.FontSize = 24;
% ax.XLim=([1 ncycles]);
% ax.YLim=([min(eff_nobact)-1  max(eff_nobact)]);
% legend({'H2 , msalt=0'},...
%     'FontSize',20,'TextColor','black','Location','west')
% legend({'pure water, no diffusion, no archae','pure water, no diffusion, archae n0=1e8',...
%     'pure water, diffusion, archae n0=1e8',' pure water, no diffusion, archae n0=1e9'},...
%     'FontSize',20,'TextColor','black','Location','west')
% 
% nT=numel(states_nobact);
% nbacteria= zeros(nT,1);
% nbacteria9= zeros(nT,1);
% nbacteria9_salt= zeros(nT,1);
% nbacteria_salt= zeros(nT,1);
% ncells=G.cells.num;
% for i = 1:nT
%     nbacteria(i)=sum(states_bact{i}.nbact);
%     nbacteria_salt(i)=sum(states_bact_salt{i}.nbact);
%     nbacteria9(i)=sum(states_bact9{i}.nbact);
%     nbacteria9_salt(i)=sum(states_bact9_salt{i}.nbact);
% end
% f31=figure('Name','nbacteria','NumberTitle','off');
% f31.Position(3:4) = [900 700];
% plot(0:nT,[ncells*nbact0;nbacteria],'m','MarkerSize',11,'LineWidth',3)
% hold on;
% plot([ncells*nbact0;nbacteria_salt],'r--','MarkerSize',11,'LineWidth',2)
% hold on;
% plot([ncells*nbact0;nbacteria9],'k','MarkerSize',11,'LineWidth',3)
% hold on;
% plot([ncells*nbact0;nbacteria9_salt],'b--','MarkerSize',11,'LineWidth',2)
% title('Total methanogenic Archae population','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'N_{archae}'},'FontWeight','bold','Color','k')
% ax = gca; 
% ax.XLim=([0 nT]);
% ax.YLim=([0  max(nbacteria)+1e3]);
% ax.FontSize = 24;
% legend({'n0=1.e8, pure water','n0=1.e8, salt water','n0=1.e9, pure water',...
%    'n0=1.e9, salt water'},...
%    'FontSize',20,'TextColor','black','Location','west')
% 

% nT=numel(states_nobact);
% f20=figure('Name','loss-gain microbial activity','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:nT,H2_loss_bact,'m','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,H2_loss_bact_salt,'r--','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:nT,H2_loss_bact9,'k','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,H2_loss_bact9_salt,'b--','MarkerSize',11,'LineWidth',2)
% title(' H2 loss due to microbial activity','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.XMinorTick='on';
% ax.FontSize = 24;
% ax.XLim=([0 nT]);
% %ax.YLim=([min(CO2_loss_bact) max(C1_gain_bact)]);
% % legend({'H2 consumption','CO2 consumption',...
% %    'C1 gain'},...
% %    'FontSize',20,'TextColor','black','Location','west')
% legend({'n0=1.e8, pure water','n0=1.e8, salt water','n0=1.e9, pure water',...
%    'n0=1.e9, salt water'},...
%    'FontSize',20,'TextColor','black','Location','west')
% 


% nT=numel(states_nobact);
% f20=figure('Name','loss-gain microbial activity','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:nT,CO2_loss_bact,'m','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,CO2_loss_bact_salt,'r--','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:nT,CO2_loss_bact9,'k','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,CO2_loss_bact9_salt,'b--','MarkerSize',11,'LineWidth',2)
% title(' CO2 loss due to microbial activity','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'CO2 consumption (%)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.XMinorTick='on';
% ax.FontSize = 24;
% ax.XLim=([0 nT]);
% %ax.YLim=([min(CO2_loss_bact) max(C1_gain_bact)]);
% % legend({'H2 consumption','CO2 consumption',...
% %    'C1 gain'},...
% %    'FontSize',20,'TextColor','black','Location','west')
% legend({'n0=1.e8, pure water','n0=1.e8, salt water','n0=1.e9, pure water',...
%    'n0=1.e9, salt water'},...
%    'FontSize',20,'TextColor','black','Location','west')
% 
% 
% 
% 
% nT=numel(states_nobact);
% f20=figure('Name','loss-gain microbial activity','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:nT,C1_gain_bact,'m','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,C1_gain_bact_salt,'r--','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:nT,C1_gain_bact9,'k','MarkerSize',11,'LineWidth',3)
% hold on;
% plot(1:nT,C1_gain_bact9_salt,'b--','MarkerSize',11,'LineWidth',2)
% title(' C1 gain due to microbial activity','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'C1 gain (%)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.XMinorTick='on';
% ax.FontSize = 24;
% ax.XLim=([0 nT]);
% %ax.YLim=([min(CO2_loss_bact) max(C1_gain_bact)]);
% % legend({'H2 consumption','CO2 consumption',...
% %    'C1 gain'},...
% %    'FontSize',20,'TextColor','black','Location','west')
% legend({'n0=1.e8, pure water','n0=1.e8, salt water','n0=1.e9, pure water',...
%    'n0=1.e9, salt water'},...
%    'FontSize',20,'TextColor','black','Location','west')


% 
nT=numel(states_nobact);
f20=figure('Name','H2totmass with microbial activity','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_loss_bact,'b','MarkerSize',11,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_bact,'r--','MarkerSize',11,'LineWidth',2)
hold on;
plot(1:nT,C1_gain_bact,'k--','MarkerSize',11,'LineWidth',2)
title(' C1gain with microbial activity','FontSize',24,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 24;
ax.XLim=([0 nT]);
%ax.YLim=([0 max(CO2_loss_bact)+1]);
%legend({'H2 , msalt=0'},...
 %   'FontSize',2,'TextColor','black','Location','west')
legend({'pure water, C1',...
   'pure water, abs(C1)'},...
   'FontSize',20,'TextColor','black','Location','west')










% 
% f20=figure('Name','H2totalMass','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:nT,results_nobact.totMassH2,'b','MarkerSize',7,'LineWidth',2)
% hold on;
% plot(1:nT,results_nobact_moldiff.totMassH2,'r--','MarkerSize',7,'LineWidth',2)
% hold on;
% plot(1:nT,results_nobact_disp.totMassH2,'k--','MarkerSize',7,'LineWidth',2)
% title(' H2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'H2 mass (kg)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.FontSize = 16;
% legend({'H2 total mass, no diff, no disp','H2 total mass, diff','H2 total mass, diff and disp'},...
%     'FontSize',16,'TextColor','black','Location','west')
% 
% 
% nT=numel(states_nobact);
% f20=figure('Name','C1_well','NumberTitle','off');
% f20.Position(3:4) = [900 700];
% plot(1:nT,results_nobact.C1_well./schedule.step.val,'b','MarkerSize',11,'LineWidth',2)
% hold on;
% plot(1:nT,results_bact.C1_well./schedule.step.val,'r--','MarkerSize',11,'LineWidth',2)
% title(' C1 well rate in pure water','FontSize',24,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'C1 Production (kg/day)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.XMinorTick='on';
% ax.YMinorTick='on';
% ax.FontSize = 24; 
% legend({'No archae','Archae'},...
%     'FontSize',20,'TextColor','black','Location','west')
% % 


% 
% 
% 
% 
% 
% %% --- Plot Benchmark Results ---
% plotBenchmarckAEGE2023(numel(states_nobact), ...
%     results_nobact.pressure, results_bact.pressure, ...
%     H2_loss, results_nobact.totMassH2, results_bact.totMassH2, ...
%     G, states_bact, nbact0);

%% Post-processing Function
function results = postProcessResults(states, ws, model, caseType)
% Post-process simulation results
%
% Parameters:
%   states - Cell array of state variables
%   ws - Well solution data
%   model - Simulation model
%   caseType - 'bact' or 'nobact'

% Get component indices
namecp = model.EOSModel.getComponentNames();
ncomp = model.EOSModel.getNumberOfComponents();
nT = numel(states);

indH2 = find(strcmp(namecp, 'H2'));
indCO2 = find(strcmp(namecp, 'CO2'));
indC1 = find(strcmp(namecp, 'C1'));

% Initialize result structure
results = struct();
resultFields = {'xH2', 'yH2', 'xCO2', 'yCO2', 'yC1', 'pressure', ...
    'H2_well', 'CO2_well', 'C1_well', 'totMassH2', ...
    'totMassCO2', 'totMassC1', 'FractionMassH2', ...
    'FractionMassCO2', 'FractionMassC1', 'totMassComp'};

for i = 1:numel(resultFields)
    results.(resultFields{i}) = zeros(nT, 1);
end

% Process results for each time step
for i = 1:nT
    results.xH2(i) = max(states{i}.x(:, indH2));
    results.yH2(i) = max(states{i}.y(:, indH2));
    results.yC1(i) = max(states{i}.y(:, indC1));
    results.xCO2(i) = max(states{i}.x(:, indCO2));
    results.yCO2(i) = max(states{i}.y(:, indCO2));
    results.pressure(i) = mean(states{i}.pressure(:));
    results.H2_well(i) = ws{i}.H2;
    results.CO2_well(i) = ws{i}.CO2;
    results.C1_well(i) = ws{i}.C1;

    % Calculate mass balances
    for j = 1:ncomp
        results.totMassComp(i) = results.totMassComp(i) + ...
            sum(states{i}.FlowProps.ComponentTotalMass{j});
    end

    results.totMassH2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});
    results.totMassCO2(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCO2});
    results.totMassC1(i) = sum(states{i}.FlowProps.ComponentTotalMass{indC1});

    results.FractionMassH2(i) = results.totMassH2(i) / results.totMassComp(i);
    results.FractionMassCO2(i) = results.totMassCO2(i) / results.totMassComp(i);
    results.FractionMassC1(i) = results.totMassC1(i) / results.totMassComp(i);
end

results.caseType = caseType;
end

%% Efficiency Calculation Function
% function efficiency = calculateH2Efficiency(H2_well, nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles)
% % Calculate H2 production efficiency
% %
% % Parameters:
% %   H2_well - H2 well data over time
% %   nbuildUp, nrest, ninject, nidle, nprod - Time step counts
% %   ncycles - Number of cycles
% 
% ndeb0 = nbuildUp + nrest;
% njcycle = ninject + nidle + nprod + nidle1;
% 
% mH2_injected = sum(H2_well(1:nbuildUp));
% mH2_produced = 0.0;
% 
% for cycle = 1:ncycles
%     ndebi = ndeb0 + (cycle-1)*njcycle;
%     nj1 = ndebi + ninject;
% 
%     mH2_injected = mH2_injected + sum(H2_well(ndebi+1:nj1));
%     mH2_produced = mH2_produced + sum(H2_well(nj1+nidle+1:nj1+nidle+nprod));
% end
% 
% efficiency = abs(mH2_produced / mH2_injected) * 100;
% end






function efficiency_H2_cycle = calculateH2Efficiency(H2_well, nbuildUp, nrest, ninject, nidle, nprod, nidle1, ncycles)
% Calculate H2 production efficiency
%
% Parameters:
%   H2_well - H2 well data over time
%   nbuildUp, nrest, ninject, nidle, nprod - Time step counts
%   ncycles - Number of cycles


ndeb0=nbuildUp+nrest;
njcycle=ninject+nidle+nprod+nidle1;
mH2_injected=sum(H2_well(1:nbuildUp));
mH2_produced=0.0;
efficiency_H2_cycle=zeros(ncycles,1);

for cycle=1:ncycles
    ndebi=ndeb0+(cycle-1)*njcycle; nj1=ndebi+ninject;
    mH2_injected=mH2_injected+sum(H2_well(ndebi+1:nj1));
    mH2_produced=mH2_produced+sum(H2_well(nj1+nidle+1:nj1+nidle+nprod));

    efficiency_H2_cycle(cycle)=abs(mH2_produced./mH2_injected).*100;
end

end












%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>