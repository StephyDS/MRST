%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of 
%a archaea.
%This test case comes from a Benchmark in EAGE 2023
%TO DO: ERROR BUG WITH BACTERIAL=FALSE AND PRODUCTION WELL. DIVERGES WITH
%BACTERIA
% Clear workspace and initialize MRST modules
clear; clc;
%mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on 
biochemistrymodel=true;
%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(31,31,8);  % Grid cells in x, y, z directions
[Lx,Ly,Lz] = deal(1525,1525,50);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, Ly, Lz];


% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 1168;                % Reservoir depth in meters
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth_res;
G = computeGeometry(G);

% Define rock properties
K=[100, 100, 10].*milli*darcy;
rock = makeRock(G, K, 0.2);  % Default permeability and porosity


%% Fluid Properties Initialization
% Define compositional fluid model (with CoolProp library support)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'C1'});

% Fluid density and viscosity (kg/m^3 and cP)
[rhow, rhog] = deal(996.52 * kilogram / meter^3, 69.974 * kilogram / meter^3);
[viscow, viscog] = deal(0.65403 * centi * poise, 0.013933 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(4.5157e-5, 1.09e-2 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.0, 0.0);
P0=82 * barsa;
T0 = 313.15;                % Initial temperature (K)
fluid = initSimpleADIFluid('phases', 'OG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', P0, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcOG = @(sw) Pe * sw.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Simulation Parameters
% Set total time, pore volume, and injection rate
niter=120;
TotalTime = niter*day;
rate = 1e6*meter^3/day; 


%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);
deltaT =rampupTimesteps(TotalTime, 1*day, 0);
schedule = simpleSchedule(deltaT);
nj1=90;nj2=100;nj3=110;
schedule.step.control(1:nj1)=1;
schedule.step.control(nj1+1:nj2)=2;
schedule.step.control(nj2+1:nj3)=3;
schedule.step.control(nj3+1:end)=4;

%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
W3 = [];
W4 = [];
tmp = cell(4,1);
n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
schedule.control = struct('W',tmp);

cellInd =[481;1442;2403;3364;4325;5286;6247;7208];


% Add a production well at the identified cells with specified properties
W0 = addWell([], G, rock, cellInd, ...
    'Name', 'Build-up', ...                       % Well name
    'Type', 'rate', ...                      % Well type (rate control)
    'Val', rate, ...           % Production rate
    'Sign', 1, ...               % Sign
    'Compi', [0, 1]);                        % Component indices
W0(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});
W0(1).T = T0;
% Injection well parameters
W1 = addWell([], G, rock, cellInd, ...
    'Name', 'Injector', ...                       % Well name
    'Type', 'rate', ...                      % Well type (rate control)
    'Sign', 1, ...               % Sign
    'Val', 0.5*rate, ...           % Production rate
    'Compi', [0, 1]);
W1(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});
W1(1).T = T0;
%Idle period
W2 = addWell([], G, rock, cellInd, ...
    'Name', 'Rest', ...                       % Well name
    'Type', 'rate', ...                      % Well type (rate control)
    'Val', 0, ...           % Production rate
    'Sign', 0, ...               % Sign
    'Compi', [0, 1]);
W2(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period
W2(1).T = T0;
%production
Pwell=40*barsa; 
W3 = addWell([], G, rock, cellInd, ...
    'Name', 'Production', ...                       % Well name
    'Type', 'bhp', ...                      % Well type (rate control)
    'Val', Pwell, ...           % Production rate
    'Sign', -1, ...               % Sign
    'Compi', [0, 1]);
W3(1).components = [0.0, 0.95,  0.05, 0.0];  %production
%W3.lims.bhp= P0;qq
W3(1).T = T0;
%Idle period
W4 = addWell([], G, rock, cellInd, ...
    'Name', 'Idle', ...                       % Well name
    'Type', 'rate', ...                      % Well type (rate control)
    'Val', 0, ...           % Production rate
    'Sign', 0, ...               % Sign
    'Compi', [0, 1]);
W4(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period
W4(1).T = T0;

% schedule.control(1).W = W1;
% schedule.control(2).W = W2;
% schedule.control(3).W = W3;
% schedule.control(4).W = W4;

schedule = createCyclicScenario(320*day, 1*day, 4, 180, 5, 15, 15, [W0;W2;W1;W3]);
%% Model Setup: Compositional Model with Bacterial Growth
if biochemistrymodel
    eosname='PR';% 'pr';
    eosmodel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
    %includeWater=true
    arg = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel', true,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    model.EOSModel.msalt=0;
else
    arg = {G, rock, fluid, compFluid, 'water', true, 'oil', false, 'gas', true, ...
      'liquidPhase', 'W', 'vaporPhase', 'G'};
    model = GenericOverallCompositionModel(arg{:});
end



%% Initial Conditions
% Temperature and initial saturations
T0 = 313.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.2-0.001-0.05, 0.001, 0.05, 0.8];  % Initial composition: H2O, H2, CO2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration

if biochemistrymodel
    if model.bacteriamodel
        nbact0 = 10^6;
        state0 = initCompositionalStateBacteria(model, Phydro0, T0, s0, ...
            z0, nbact0,eosmodel);
    else
        state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
    end

else
    state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
end


%% Run simulation
%mrstModule add mpfa
%model_mpfa = setMPFADiscretization(model);
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);
mrstModule add mpfa
model_mpfa = setMPFADiscretization(model);
[wellSols,states,report]= simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);


%% Plotting Results
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
nT=numel(states);
xH2=zeros(nT,1);
yH2=zeros(nT,1);
yCO2= zeros(nT,1);
for i = 1:nT
    xH2(i)=max(states{i}.x(:,indH2));
    yH2(i)=max(states{i}.y(:,indH2));
    yCO2(i)=max(states{i}.y(:,indCO2));
end

for i = 1:nT
    figure(1); clf; 
    plot(1:nT,yH2,'b')
    hold on
    plot(1:nT,yCO2,'k-')
end
title('Molar fractions')
xlabel('Time (days)')
ylabel('molar fraction')
legend('yH2','yCO2')

figure(2),clf
for i=1:nT
    clf;
    plotCellData(G,states{i}.nbact./nbact0);
    colorbar; 
    axis equal
    axis ([0 Lx  0 Ly depth_res depth_res+Lz])
    view(0,-90)
    pause(0.1)   
end
title('Archea density')
