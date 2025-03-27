%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','CH4') and The microbial activity of 
%a archaea.
%This test case comes from a Benchmark in EAGE 2023
% Clear workspace and initialize MRST modules
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on 
biochemistrymodel=true;%false; 

%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(31,31,8);  % Grid cells in x, y, z directions
[Lx,Ly,Lz] = deal(1525,1525,50);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, Ly, Lz];


% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 1000;   % Reservoir depth in meters
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
[rhow, rhog] = deal(999.7 * kilogram / meter^3, 1.2243 * kilogram / meter^3);
[viscow, viscog] = deal(1.3059 * centi * poise, 0.01763 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(5.0015e-5/ barsa, 1.0009 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.2, 0.05);
P0=100 * barsa;
fluid = initSimpleADIFluid('phases', 'OG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', P0, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcOG = @(sw) Pe * sw.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Simulation Parameters
% Set total time, pore volume, and injection rate
niter=170; 
TotalTime = niter*day;
rate = 1e6*meter^3/day; 


%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);
deltaT =rampupTimesteps(TotalTime, 1*day, 0);
schedule = simpleSchedule(deltaT);
nj1=60;nj2=80;nj3=110;nj4=130;nj5=160;
schedule.step.control(1:nj1)=1;
schedule.step.control(nj1+1:nj2)=2;
schedule.step.control(nj2+1:nj3)=3;
schedule.step.control(nj3+1:nj4)=4;
schedule.step.control(nj4+1:nj5)=5;
schedule.step.control(nj5+1:end)=6;

%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
W3 = [];
W4 = [];
W5 = [];
W6 = [];

tmp = cell(10,1);
n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
schedule.control = struct('W',tmp);
cellInd=zeros(nz-1,1);
for k=2:nz
    cellInd(k-1)=(k-1)*nx*ny+(n2-1)*nx+n1;
end
% Build-in well parameters
W1 = verticalWell(W1, G, rock, n1, n2, 1:nz-1, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W1(1).components = [0.0, 0.6,  0.4, 0.0];  % CO2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});

%Rest period
W2 = verticalWell(W2, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W2(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period


% Injection well parameters
W3 = verticalWell(W3, G, rock, n1, n2, 1:nz-1, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W3(1).components = [0.0, 0.95,  0.05, 0.0];  % H2-rich injection   {'H2O', 'H2', 'CO2', 'C1'});


%Idle period
W4 = verticalWell(W4, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W4(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period

%production
%Pwell=66*barsa;
%W3 = verticalWell(W3, G, rock, n1, n2, 1:nz-1, 'comp_i', [0, 1], 'Radius', 0.5, ...
 %                 'name', 'Prod', 'type', 'bhp', 'Val', Pwell, 'sign', -1);
W5 = verticalWell(W5, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                  'name', 'Prod', 'type', 'rate', 'Val', -rate, 'sign', -1);
W5(1).components = [0.0, 0.95,  0.05, 0.0];  %production

%Idle period
W6 = verticalWell(W6, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                 'name', 'Rest', 'type', 'rate', 'Val', 0.0, 'sign', 1);
W6(1).components = [0.0, 0.95,  0.05, 0.0];  % rest period


schedule.control(1).W = W1;
schedule.control(2).W = W2;
schedule.control(3).W = W3;
schedule.control(4).W = W4;
schedule.control(5).W = W5;
schedule.control(6).W = W6;

%% Model Setup: Compositional Model with Bacterial Growth
if biochemistrymodel
    eosname='sw';
    model.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);

    arg = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel', true,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    model.EOSModel.msalt=3.0; %0;
else
    eosname='pr';
    arg = {G, rock, fluid, compFluid, 'water', false, 'oil', true, 'gas', false, ...
      'liquidPhase', 'O', 'vaporPhase', 'G'};
    model = GenericOverallCompositionModel(arg{:});
end



%% Initial Conditions
% Temperature and initial saturations
T0 = 40+273.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.7, 0.0, 0.02, 0.28];  % Initial composition: H2O, H2, CO2, CH4
%Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration

if biochemistrymodel
    if model.bacteriamodel
        nbact0 = 50; 
        state0 = initCompositionalStateBacteria(model, P0, T0, s0, ...
            z0, nbact0,model.EOSModel);
    else
        state0 = initCompositionalState(model, P0, T0, s0, z0);
    end

else
    state0 = initCompositionalState(model, P0, T0, s0, z0);
end


%% ===========Compar simulation results=============================
%% open directories%%%%%%%%%%%%%%%%%%
    dir='/home/sdelage2/PROJETS/gdr_h2/MRST2024/MRST/output';
    handler_nbs0 = ResultHandler('dataDirectory',dir,'dataFolder','Benchmark2023AEGE_nbs0_170');
    handler_nbs3 = ResultHandler('dataDirectory',dir,'dataFolder','Benchmark2023AEGE_nbs3_170');
    handler_bs0 = ResultHandler('dataDirectory',dir,'dataFolder','Benchmark2023AEGE_bs0_170');
    handler_bs3 = ResultHandler('dataDirectory',dir,'dataFolder','Benchmark2023AEGE_bs3_170');
    m = handler_nbs0.numelData();
    states_nbs0 = cell(m, 1);
    states_nbs3 = cell(m, 1);
    states_bs0 = cell(m, 1);
    states_bs3 = cell(m, 1);
    for i = 1:m
        states_nbs0{i} = handler_nbs0{i};
        states_nbs3{i} = handler_nbs3{i};
        states_bs0{i} = handler_bs0{i};
        states_bs3{i} = handler_bs3{i};
    end

    %% initialise unknowns %%%%%%%%%%%%%%%%%%
    namecp = model.EOSModel.getComponentNames();
    ncomp=model.EOSModel.getNumberOfComponents;

    indH2=find(strcmp(namecp,'H2'));
    indCO2= find(strcmp(namecp,'CO2'));
    indCH4= find(strcmp(namecp,'C1'));
    nT=numel(states_nbs0);

    xH2_nbs0=zeros(nT,1);
    yH2_nbs0=zeros(nT,1);
    xCO2_nbs0= zeros(nT,1);
    yCO2_nbs0= zeros(nT,1);
    yCH4_nbs0= zeros(nT,1);
    pressure_nbs0=zeros(nT,1);
    H2_well_nbs0=zeros(nT,1);
    H2_well_bs0=zeros(nT,1);
    CO2_well_nbs0=zeros(nT,1);
    CO2_well_bs0=zeros(nT,1);
    CH4_well_nbs0=zeros(nT,1);
    CH4_well_bs0=zeros(nT,1);

    xH2_nbs3=zeros(nT,1);
    yH2_nbs3=zeros(nT,1);
    xCO2_nbs3= zeros(nT,1);
    yCO2_nbs3= zeros(nT,1);
    yCH4_nbs3= zeros(nT,1);
    pressure_nbs3=zeros(nT,1);
    H2_well_nbs3=zeros(nT,1);
    H2_well_bs3=zeros(nT,1);
    CO2_well_nbs3=zeros(nT,1);
    CO2_well_bs3=zeros(nT,1);
    CH4_well_nbs3=zeros(nT,1);
    CH4_well_bs3=zeros(nT,1);

    xH2_bs0=zeros(nT,1);
    yH2_bs0=zeros(nT,1);
    xCO2_bs0= zeros(nT,1);
    yCO2_bs0= zeros(nT,1);
    yCH4_bs0= zeros(nT,1);
    pressure_bs0=zeros(nT,1);

    xH2_bs3=zeros(nT,1);
    yH2_bs3=zeros(nT,1);
    xCO2_bs3= zeros(nT,1);
    yCO2_bs3= zeros(nT,1);
    yCH4_bs3= zeros(nT,1);
    pressure_bs3=zeros(nT,1);

    totMassH2_nbs0= zeros(nT,1);
    totMassCO2_nbs0= zeros(nT,1);
    totMassCH4_nbs0= zeros(nT,1);
    FractionMassCO2_nbs0= zeros(nT,1);
    FractionMassH2_nbs0= zeros(nT,1);
    FractionMassCH4_nbs0= zeros(nT,1);
    
    totMassComp_nbs0= zeros(nT,1);
    totMassH2_nbs3= zeros(nT,1);
    totMassCO2_nbs3= zeros(nT,1);
    totMassCH4_nbs3= zeros(nT,1);
    FractionMassCO2_nbs3= zeros(nT,1);
    FractionMassH2_nbs3= zeros(nT,1);
    FractionMassCH4_nbs3= zeros(nT,1);
    totMassComp_nbs3= zeros(nT,1);

    totMassH2_bs0= zeros(nT,1);
    totMassCO2_bs0= zeros(nT,1);
    totMassCH4_bs0= zeros(nT,1);
    FractionMassCO2_bs0= zeros(nT,1);
    FractionMassH2_bs0= zeros(nT,1);
    FractionMassCH4_bs0= zeros(nT,1);
    totMassComp_bs0= zeros(nT,1);

    totMassH2_bs3= zeros(nT,1);
    totMassCO2_bs3= zeros(nT,1);
    totMassCH4_bs3= zeros(nT,1);
    FractionMassCO2_bs3= zeros(nT,1);
    FractionMassH2_bs3= zeros(nT,1);
    FractionMassCH4_bs3= zeros(nT,1);
    totMassComp_bs3= zeros(nT,1);


    for i = 1:nT
        xH2_nbs0(i)=max(states_nbs0{i}.x(:,indH2));
        yH2_nbs0(i)=max(states_nbs0{i}.y(:,indH2));
        yCH4_nbs0(i)=max(states_nbs0{i}.y(:,indCH4));
        xCO2_nbs0(i)=max(states_nbs0{i}.x(:,indCO2));
        yCO2_nbs0(i)=max(states_nbs0{i}.y(:,indCO2));
        pressure_nbs0(i)=mean(states_nbs0{i}.pressure(:));
        H2_well_nbs0(i)=states_nbs0{i}.wellSol.H2;
        H2_well_bs0(i)=states_bs0{i}.wellSol.H2;
        CO2_well_nbs0(i)=states_nbs0{i}.wellSol.CO2;
        CO2_well_bs0(i)=states_bs0{i}.wellSol.CO2;
        CH4_well_nbs0(i)=states_nbs0{i}.wellSol.C1;
        CH4_well_bs0(i)=states_bs0{i}.wellSol.C1;

        xH2_nbs3(i)=max(states_nbs3{i}.x(:,indH2));
        yH2_nbs3(i)=max(states_nbs3{i}.y(:,indH2));
        yCH4_nbs3(i)=max(states_nbs3{i}.y(:,indCH4));
        xCO2_nbs3(i)=max(states_nbs3{i}.x(:,indCO2));
        yCO2_nbs3(i)=max(states_nbs3{i}.y(:,indCO2));
        pressure_nbs3(i)=mean(states_nbs3{i}.pressure(:));
         H2_well_nbs3(i)=states_nbs3{i}.wellSol.H2;
        H2_well_bs3(i)=states_bs3{i}.wellSol.H2;
        CO2_well_nbs3(i)=states_nbs3{i}.wellSol.CO2;
        CO2_well_bs3(i)=states_bs3{i}.wellSol.CO2;
        CH4_well_nbs3(i)=states_nbs3{i}.wellSol.C1;
        CH4_well_bs3(i)=states_bs3{i}.wellSol.C1;

        xH2_bs0(i)=max(states_bs0{i}.x(:,indH2));
        yH2_bs0(i)=max(states_bs0{i}.y(:,indH2));
        yCH4_bs0(i)=max(states_bs0{i}.y(:,indCH4));
        xCO2_bs0(i)=max(states_bs0{i}.x(:,indCO2));
        yCO2_bs0(i)=max(states_bs0{i}.y(:,indCO2));
        pressure_bs0(i)=mean(states_bs0{i}.pressure(:));

        xH2_bs3(i)=max(states_bs3{i}.x(:,indH2));
        yH2_bs3(i)=max(states_bs3{i}.y(:,indH2));
        yCH4_bs3(i)=max(states_bs3{i}.y(:,indCH4));
        xCO2_bs3(i)=max(states_bs3{i}.x(:,indCO2));
        yCO2_bs3(i)=max(states_bs3{i}.y(:,indCO2));
        pressure_bs3(i)=mean(states_bs3{i}.pressure(:));
    end

    %% mass calculations
    for i = 1:nT
        for j=1:ncomp
            totMassComp_nbs0(i)=totMassComp_nbs0(i)+sum(states_nbs0{i}.FlowProps.ComponentTotalMass{j});
            totMassComp_nbs3(i)=totMassComp_nbs3(i)+sum(states_nbs3{i}.FlowProps.ComponentTotalMass{j});
            totMassComp_bs0(i)=totMassComp_bs0(i)+sum(states_bs0{i}.FlowProps.ComponentTotalMass{j});
            totMassComp_bs3(i)=totMassComp_bs3(i)+sum(states_bs3{i}.FlowProps.ComponentTotalMass{j});
        end
        totMassH2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_nbs0(i)=totMassH2_nbs0(i)/totMassComp_nbs0(i);
        FractionMassCO2_nbs0(i)=totMassCO2_nbs0(i)/totMassComp_nbs0(i);
        FractionMassCH4_nbs0(i)=totMassCH4_nbs0(i)/totMassComp_nbs0(i);

        totMassH2_nbs3(i)=sum(states_nbs3{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_nbs3(i)=sum(states_nbs3{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_nbs3(i)=sum(states_nbs3{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_nbs3(i)=totMassH2_nbs3(i)/totMassComp_nbs3(i);
        FractionMassCO2_nbs3(i)=totMassCO2_nbs3(i)/totMassComp_nbs3(i);
        FractionMassCH4_nbs3(i)=totMassCH4_nbs3(i)/totMassComp_nbs3(i);

        totMassH2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_bs0(i)=totMassH2_bs0(i)/totMassComp_bs0(i);
        FractionMassCO2_bs0(i)=totMassCO2_bs0(i)/totMassComp_bs0(i);
        FractionMassCH4_bs0(i)=totMassCH4_bs0(i)/totMassComp_bs0(i);

        totMassH2_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_bs3(i)=totMassH2_bs3(i)/totMassComp_bs3(i);
        FractionMassCO2_bs3(i)=totMassCO2_bs3(i)/totMassComp_bs3(i);
        FractionMassCH4_bs3(i)=totMassCH4_bs3(i)/totMassComp_bs3(i);  
    end

    %% Calculate H2 production efficiency 
    %% in the reservoir
    mH2_injected_nbs0=totMassH2_nbs0(nj1)-totMassH2_nbs0(1)...
        +totMassH2_nbs0(nj3)-totMassH2_nbs0(nj2+1);
    mH2_produced_nbs0=totMassH2_nbs0(nj4+1)-totMassH2_nbs0(nj5);
    mH2_injected_nbs3=totMassH2_nbs3(nj1)-totMassH2_nbs3(1)...
        +totMassH2_nbs3(nj3)-totMassH2_nbs3(nj2+1);
    mH2_produced_nbs3=totMassH2_nbs3(nj4+1)-totMassH2_nbs3(nj5);
    mH2_injected_bs0=totMassH2_bs0(nj1)-totMassH2_bs0(1)...
        +totMassH2_bs0(nj3)-totMassH2_bs0(nj2+1);
    mH2_produced_bs0=totMassH2_bs0(nj4+1)-totMassH2_bs0(nj5);
    mH2_injected_bs3=totMassH2_bs3(nj1)-totMassH2_bs3(1)...
        +totMassH2_bs3(nj3)-totMassH2_bs3(nj2+1);
    mH2_produced_bs3=totMassH2_bs3(nj4+1)-totMassH2_bs3(nj5);
   
    Efficiency_H2_nbs0=(mH2_produced_nbs0/mH2_injected_nbs0) * 100;
    Efficiency_H2_nbs3=(mH2_produced_nbs3/mH2_injected_nbs3) * 100;
    Efficiency_H2_bs0=(mH2_produced_bs0/mH2_injected_bs0) * 100;
    Efficiency_H2_bs3=(mH2_produced_bs3/mH2_injected_bs3) * 100;

    %% Display H2 production efficiency
    fprintf('H2 Production Efficiency_nobact, msalt=0 : %.2f%%\n', Efficiency_H2_nbs0);
    fprintf('H2 Production Efficiency_nobact, msalt=3 : %.2f%%\n', Efficiency_H2_nbs3);
    fprintf('H2 Production Efficiency_bact, msalt=0 : %.2f%%\n', Efficiency_H2_bs0);
    fprintf('H2 Production Efficiency_bact, msalt=3 : %.2f%%\n', Efficiency_H2_bs3);

    %% in the wellmH2_injected_bs0=sum(H2_well_bs0(1:60))+sum(H2_well_bs0(81:110))
    mH2_well_injected_nbs0=sum(H2_well_nbs0(1:nj1))+sum(H2_well_nbs0(nj2+1:nj3));
    mH2_well_injected_bs0=sum(H2_well_bs0(1:nj1))+sum(H2_well_bs0(nj2+1:nj3));
    mH2_well_injected_nbs3=sum(H2_well_nbs3(1:nj1))+sum(H2_well_nbs3(nj2+1:nj3));
    mH2_well_injected_bs3=sum(H2_well_bs3(1:nj1))+sum(H2_well_bs3(nj2+1:nj3));
    mH2_well_produced_nbs0=sum(H2_well_nbs0(nj4+1:nj5));
    mH2_well_produced_bs0=sum(H2_well_bs0(nj4+1:nj5));
    mH2_well_produced_nbs3=sum(H2_well_nbs3(nj4+1:nj5));
    mH2_well_produced_bs3=sum(H2_well_bs3(nj4+1:nj5));

    Efficiency_H2_well_nbs0=abs(mH2_well_produced_nbs0./mH2_well_injected_nbs0).*100;
    Efficiency_H2_well_bs0=abs(mH2_well_produced_bs0./mH2_well_injected_bs0).*100;
    Efficiency_H2_well_nbs3=abs(mH2_well_produced_nbs3./mH2_well_injected_nbs3).*100;
    Efficiency_H2_well_bs3=abs(mH2_well_produced_bs3./mH2_well_injected_bs3).*100;

    %% Display H2 production efficiency
    fprintf('In the well, H2 Production Efficiency_nobact, msalt=0 : %.2f%%\n', Efficiency_H2_well_nbs0);
    fprintf('In the well, H2 Production Efficiency_nobact, msalt=3 : %.2f%%\n', Efficiency_H2_well_nbs3);
    fprintf('In the well, H2 Production Efficiency_bact, msalt=0 : %.2f%%\n', Efficiency_H2_well_bs0);
    fprintf('In the well, H2 Production Efficiency_bact, msalt=3 : %.2f%%\n', Efficiency_H2_well_bs3);
   
    %% Calculate percentage of H2 loss
    H2_loss_percentage_nbs0 = (abs(totMassH2_nbs0-totMassH2_bs0)./totMassH2_nbs0) * 100;
    H2_loss_percentage_nbs3 = (abs(totMassH2_nbs3-totMassH2_bs3)./totMassH2_nbs3) * 100;
    %% Calculate percentage of CO2 loss
    CO2_loss_percentage_nbs0 = (abs(totMassCO2_nbs0-totMassCO2_bs0)./totMassCO2_nbs0) * 100;
    CO2_loss_percentage_nbs3 = (abs(totMassCO2_nbs3-totMassCO2_bs3)./totMassCO2_nbs3) * 100;
    %% Calculate percentage of CH4 production
    CH4_loss_percentage_nbs0 = (abs(totMassCH4_nbs0-totMassCH4_bs0)./totMassCH4_nbs0) * 100;
    CH4_loss_percentage_nbs3 = (abs(totMassCH4_nbs3-totMassCH4_bs3)./totMassCH4_nbs3) * 100;

%% Display final H2 loss
fprintf('Total H2 loss due to bacterial effects, msalt=0: %.2f%%\n', H2_loss_percentage_nbs0(end));
fprintf('Total CO2 loss due to bacterial effects, msalt=0: %.2f%%\n', CO2_loss_percentage_nbs0(end));
fprintf('Total CH4 production due to bacterial effects, msalt=0: %.2f%%\n', CH4_loss_percentage_nbs0(end));
fprintf('Total H2 loss due to bacterial effects, msalt=3: %.2f%%\n', H2_loss_percentage_nbs3(end));
fprintf('Total CO2 loss due to bacterial effects, msalt=3: %.2f%%\n', CO2_loss_percentage_nbs3(end));
fprintf('Total CH4 production due to bacterial effects, msalt=3: %.2f%%\n', CH4_loss_percentage_nbs3(end));
 
%% plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
%% plot pressure
f01=figure('Name','Pressure_compar_bact_msalt0','NumberTitle','off');
f01.Position(3:4) = [900 700];
plot(1:nT,pressure_bs0./1e5,'b','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_nbs0./1e5,'r--','MarkerSize',8,'LineWidth',2)
title('Mean pressure in the pure water reservoir','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'pressure (bar)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'pressure, with archae','pressure, no archae'},...
    'FontSize',16,'TextColor','black',...
    'Location','best')

f02=figure('Name','Pressure_compar_bact_msalt3','NumberTitle','off');
f02.Position(3:4) = [900 700];
plot(1:nT,pressure_bs3./1e5,'b','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_nbs3./1e5,'r--','MarkerSize',8,'LineWidth',2)
title('Mean pressure in the salt water reservoir','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'pressure (bar)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'pressure, with archae','pressure, no archae'},...
    'FontSize',16,'TextColor','black',...
    'Location','best')

%% plot molar H2 fraction in Gaz
f11=figure('Name','maxyH2_compar_bact_msalt0','NumberTitle','off');
f11.Position(3:4) = [900 700];
plot(0:nT,[0;yH2_bs0],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yH2_nbs0],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum H_2 Molar fractions in gaz, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{g,H_2}, with archae','X_{g,H_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f12=figure('Name','maxyH2_compar_bact_msalt3','NumberTitle','off');
f12.Position(3:4) = [900 700];
plot(0:nT,[0;yH2_bs3],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yH2_nbs3],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum H_2 Molar fractions in gaz, salt water','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{g,H_2}, with archae','X_{g,H_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

%% plot H2 molar fraction in Liquid
f13=figure('Name','maxxH2_compar_bact_msalt0','NumberTitle','off');
f13.Position(3:4) = [900 700];
plot(0:nT,[0;xH2_bs0],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xH2_nbs0],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum H_2 Molar fractions in liquid, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
ax.XTick
legend({'X_{l,H_2}, with archae','X_{l,H_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f14=figure('Name','maxxH2_compar_bact_msalt3','NumberTitle','off');
f14.Position(3:4) = [900 700];
plot(0:nT,[0;xH2_bs3],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xH2_nbs3],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum H_2 Molar fractions in liquid, salt water','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{l,H_2}, with archae','X_{l,H_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')


%% plot molar CO2 fraction in Gaz
f15=figure('Name','maxyCO2_compar_bact_msalt0','NumberTitle','off');
f15.Position(3:4) = [900 700];
plot(0:nT,[0;yCO2_bs0],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yCO2_nbs0],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum CO_2 Molar fractions in gaz, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{g,CO_2}, with archae','X_{g,CO_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f16=figure('Name','maxyCO2_compar_bact_msalt3','NumberTitle','off');
f16.Position(3:4) = [900 700];
plot(0:nT,[0;yCO2_bs3],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;yCO2_nbs3],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum CO_2 Molar fractions in gaz, salt water','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{g,CO_2}, with archae','X_{g,CO_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

%% plot CO2 molar fraction in Liquid
f17=figure('Name','maxxCO2_compar_bact_msalt0','NumberTitle','off');
f17.Position(3:4) = [900 700];
plot(0:nT,[0;xCO2_bs0],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xCO2_nbs0],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum CO_2 Molar fractions in liquid, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'CO_{l,H_2}, with archae','CO_{l,H_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f18=figure('Name','maxxCO2_compar_bact_msalt3','NumberTitle','off');
f18.Position(3:4) = [900 700];
plot(0:nT,[0;xCO2_bs3],'b','MarkerSize',7,'LineWidth',2)
hold on
plot(0:nT,[0;xCO2_nbs3],'r--','MarkerSize',8,'LineWidth',2)
title('Maximum CO_2 Molar fractions in liquid, salt water','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO_2 molar fraction '},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'X_{l,CO_2}, with archae','X_{l,CO_2}, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

%% plot well%%%%%%%%%%%%%%%%
f20=figure('Name','H2_well','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_nbs0,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on'
ax.YMinorTick='on'
ax.FontSize = 16; 
legend({'H2, with archae','H2, no archae'},...
    'FontSize',16,'TextColor','black','Location','west')

f21=figure('Name','CO2_well','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(1:nT,CO2_well_bs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_well_nbs0,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on'
ax.YMinorTick='on'
ax.FontSize = 16; 
legend({'CO2, with archae','CO2, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f22=figure('Name','C1_well','NumberTitle','off');
f22.Position(3:4) = [900 700];
plot(1:nT,CH4_well_bs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CH4_well_nbs0,'r--','MarkerSize',7,'LineWidth',2)
title(' C1 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'C1 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on'
ax.YMinorTick='on'
ax.FontSize = 16; 
legend({'C1, with archae','C1, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f23=figure('Name','H2_CO2_CH4_well','NumberTitle','off');
f23.Position(3:4) = [900 700];
plot(0:nT,[0;H2_well_nbs0],'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(0:nT,[0;CO2_well_nbs0],'r--','MarkerSize',7,'LineWidth',2)
hold on;
plot(0:nT,[0;CH4_well_nbs0],'k:','MarkerSize',7,'LineWidth',2)
title(' mass of H2, CO2 and C1 in the well','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2, CO2, C1'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'H2, no archae','CO2, no archae','C1, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f24=figure('Name','H2_well_salt','NumberTitle','off');
f24.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs3,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate, brine','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on'
ax.YMinorTick='on'
ax.FontSize = 16; 
legend({'H2, with archae','H2, no archae'},...
    'FontSize',16,'TextColor','black','Location','west')


f25=figure('Name','H2_well_bact_comparsalt','NumberTitle','off');
f25.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_bs3,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate with archae','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on'
ax.YMinorTick='on'
ax.FontSize = 16; 
legend({'H2, msalt=0','H2, msalt=3'},...
    'FontSize',16,'TextColor','black','Location','west')




% f4=figure('Name','fractionmass_compar_bact','NumberTitle','off');
% f4.Position(3:4) = [900 700];
% plot(0:nT,[0;FractionMassH2],'b','MarkerSize',7,'LineWidth',2)
% hold on
% plot(0:nT,[0;FractionMassH2_nobact],'r--','MarkerSize',8,'LineWidth',2)
% 
% title('Total H_2 mass fraction, m_{salt}=3 ','FontSize',16,'FontWeight','bold','Color','k')
% xlabel({'time (days)'},'FontWeight','bold','Color','k')
% ylabel({'H_2 mass fraction'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.FontSize = 16; 
% 
% legend({'X_{H_2}, with archae','X_{H_2}, no archae'},...
%     'FontSize',16,'TextColor','black',...
%     'Location','best')


if biochemistrymodel && model.bacteriamodel
    nbacteria_bs0= zeros(nT,1);
    nbacteria_bs3= zeros(nT,1);
    pv=model.operators.pv;
    ncells=G.cells.num;
    for i = 1:nT
       Swi_bs0 = states_bs0{i}.s(:,1);
       Swi_bs3 = states_bs3{i}.s(:,1);
       Swpvi_bs0=Swi_bs0.*pv;
       Swpvi_bs3=Swi_bs3.*pv;
       Swpv_bs0=sum(Swpvi_bs0);
       Swpv_bs3=sum(Swpvi_bs3);
       nbacteria_bs0(i)=sum(states_bs0{i}.nbact.*Swpvi_bs0)/Swpv_bs0;
       nbacteria_bs3(i)=sum(states_bs3{i}.nbact.*Swpvi_bs3)/Swpv_bs3;
    end

   f31=figure('Name','nbacteria','NumberTitle','off');
   f31.Position(3:4) = [900 700];
   plot(0:nT,[nbact0;nbacteria_bs0],'b','MarkerSize',7,'LineWidth',2)
   hold on
   plot(0:nT,[nbact0;nbacteria_bs3],'r--','MarkerSize',7,'LineWidth',2)
   title('Methanogenic Archae in pure water and brine','FontSize',16,'FontWeight','bold','Color','k')
   xlabel({'time (days)'},'FontWeight','bold','Color','k')
   ylabel({'N_{archae}'},'FontWeight','bold','Color','k')
   ax = gca;
   ax.FontSize = 16; 
   legend({'N_{archae}, pure water','N_{archae}, brine'},'FontSize',16,'TextColor','black',...
    'Location','best')

   diff_nbacteria=abs(nbacteria_bs0-nbacteria_bs3)./nbacteria_bs0;
   f32=figure('Name','diff_nbacteria','NumberTitle','off');
   f32.Position(3:4) = [900 700];
   plot(0:nT,[0;diff_nbacteria],'b','MarkerSize',7,'LineWidth',2)
   title('Impact of salinity in Archae population','FontSize',16,'FontWeight','bold','Color','k')
   xlabel({'time (days)'},'FontWeight','bold','Color','k')
   ylabel({'Archae deviation'},'FontWeight','bold','Color','k')
   ax = gca;
   ax.FontSize = 16; 
   legend({'N_{archae} deviation'},'FontSize',16,'TextColor','black',...
    'Location','best')





% f11=figure('Name','nbacteria_t5','NumberTitle','off');
% f11.Position(3:4) = [900 700];
% plotCellData(G,states{5}.nbact);
% colorbar; 
% axis equal
% axis ([0 Lx  0 Ly depth_res depth_res+Lz])
% view(0,-90)
% title('Methanogenic Archae density, 5 days','FontSize',16,'FontWeight','bold','Color','k')
% xlabel({'x (m)'},'FontWeight','bold','Color','k')
% ylabel({'y (m)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.FontSize = 16; 
% %legend({'n_{archae}'},...
%  %   'FontSize',13,'TextColor','black',...
% %    'Location','east')
% clim([0 nbact0])
end