%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of 
%a archaea.
%This test case comes from a Benchmark in EAGE 2023
% Clear workspace and initialize MRST modules
clear; clc;
%mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on 
biochemistrymodel=true;%false; 
writedatafile=true;
%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(31,31,8);  % Grid cells in x, y, z directions
[Lx,Ly,Lz] = deal(1525,1525,50);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, Ly, Lz];


% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 1000; %1310;                % Reservoir depth in meters
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
%[rhow, rhog] = deal(996.52 * kilogram / meter^3, 69.974 * kilogram / meter^3);
%[viscow, viscog] = deal(0.65403 * centi * poise, 0.013933 * centi * poise);
[rhow, rhog] = deal(999.7 * kilogram / meter^3, 1.2243 * kilogram / meter^3);
[viscow, viscog] = deal(1.3059 * centi * poise, 0.01763 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(5.0015e-5, 1.0009 / barsa);
%[cfw, cfg] = deal(4.5157e-5, 1.09e-2 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.0, 0.0);
P0=106 * barsa;%132 * barsa;
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
rate = 8e5*meter^3/day; 


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

% Injection well parameters
W1 = verticalWell(W1, G, rock, n1, n2, 1:nz, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W1(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});

%Idle period
W2 = verticalWell(W2, G, rock, n1, n2, 1:nz, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W2(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period

%production
Pwell=66*barsa;% 92*barsa; 
W3 = verticalWell(W3, G, rock, n1, n2, 1:nz, 'compi', [0, 1], 'Radius', 0.5, ...
                  'name', 'Prod', 'type', 'bhp', 'Val', Pwell, 'sign', -1);
W3(1).components = [0.0, 0.95,  0.05, 0.0];  %production

%Idle period
W4 = verticalWell(W4, G, rock, n1, n2, 1:nz, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W4(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period


schedule.control(1).W = W1;
schedule.control(2).W = W2;
schedule.control(3).W = W3;
schedule.control(4).W = W4;

%schedule = createCyclicScenario(320*day, 1*day, 4, 180, 5, 15, 15, [W0;W2;W1;W3]);
%% Model Setup: Compositional Model with Bacterial Growth
if biochemistrymodel
    eosname='sw';% 'sw';
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
    eosname='pr';
    arg = {G, rock, fluid, compFluid, 'water', false, 'oil', true, 'gas', true, ...
      'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = GenericOverallCompositionModel(arg{:});
end



%% Initial Conditions
% Temperature and initial saturations
zH2_init=0.0; %0.001;
zCO2_init=0.0; %0.001;
T0 = 313.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.2-zH2_init-zCO2_init, zH2_init, zCO2_init, 0.8];  % Initial composition: H2O, H2, CO2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration

if biochemistrymodel
    if model.bacteriamodel
        nbact0 = 1e6;  
        %model.Y_H2 = 1.7e12;  % Conversion factor for hydrogen consumption (moles/volume)
        %model.alphaH2 = 1.1e-7;
        %model.alphaCO2 = 3.2e-6;
        %model.Psigrowthmax = 1.7e-4;
        %model.b_bact = 2.3e-5/nbact0;
        state0 = initCompositionalStateBacteria(model, Phydro0, T0, s0, ...
            z0, nbact0,eosmodel);
    else
        state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
    end

else
    state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
end


%% Run simulation
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);
%mrstModule add mpfa
%model_mpfa = setMPFADiscretization(model);
%[wellSols,states,report]= simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);


%% Plotting Results
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
nT=numel(states);
xH2=zeros(nT,1);
yH2=zeros(nT,1);
xCO2= zeros(nT,1);
yCO2= zeros(nT,1);
for i = 1:nT
    xH2(i)=max(states{i}.x(:,indH2));
    yH2(i)=max(states{i}.y(:,indH2));
    xCO2(i)=max(states{i}.x(:,indCO2));
    yCO2(i)=max(states{i}.y(:,indCO2));
end

for i = 1:nT
    figure(1); clf; 
    plot(1:nT,yH2,'b')
    hold on
    plot(1:nT,yCO2,'k-')
end
title('Molar fractions in gas')
xlabel('Time (days)')
ylabel('molar fraction')
legend('yH2','yCO2')

for i = 1:nT
    figure(2); clf; 
    plot(1:nT,xH2,'b')
    hold on
    plot(1:nT,xCO2,'k-')
end
title('Molar fractions in Liquid')
xlabel('Time (days)')
ylabel('molar fraction')
legend('xH2','xCO2')

if biochemistrymodel && model.bacteriamodel
    nbacteria= zeros(nT,1);
    pv=model.operators.pv;
    ncells=G.cells.num;
    for i = 1:nT
       Swi = states{i}.s(:,1);
       Swpvi=Swi.*pv;
       Swpv=sum(Swpvi);
       nbacteria(i)=sum(states{i}.nbact.*Swpvi)/Swpv;
       %nbacteria(i)=max(states{i}.nbact);
    end

    for i = 1:nT
        figure(3); clf; 
        plot(1:nT,nbacteria,'b')
    end
    figure(4),clf
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
end

if writedatafile
    outputFileName = sprintf('H2_CO2_solubility_%s_withnbact.dat',eosname);
    fileID = fopen(outputFileName, 'wt');

    fprintf(fileID, 'iter H2_molar_fraction  CO2_molar_fraction\n');
    for i = 1:nT
       fprintf(fileID, '%d  %12.8f  %12.8f %12.8f  %12.8f\n', i, xH2(i),xCO2(i), yH2(i),yCO2(i));
    end
    fclose(fileID);

    disp(['Results written to file: ', outputFileName]);
end