%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a 3D porous medium,
% incorporating compositional fluid properties, bacterial mono modal.
% We consider a liquid phase (W) and a gas (G) phase, 4 components 
% ('H2O','H2','CO2','C1') and The microbial activity of 
%a archaea.
%This test case comes from a Benchmark in EAGE 2023
% Clear workspace and initialize MRST modules
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on 
nobact=false;%true;%

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
rate = 1e6*meter^3/day; 


%% Wells and Boundary Conditions
% Initialize wells
W1 = [];
W2 = [];
W3 = [];
W5 = [];

n1=floor(0.5*nx)+1; n2=floor(0.5*nx)+1;
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


%production
W5 = verticalWell(W5, G, rock, n1, n2, 1:nz-1, 'compi', [0, 1], 'Radius', 0.5, ...
                  'name', 'Prod', 'type', 'rate', 'Val', -rate, 'sign', -1);
W5(1).components = [0.0, 0.95,  0.05, 0.0];  %production

%% Time Stepping and Schedule
% Define schedule and solver
nls = NonLinearSolver('useRelaxation', true);

ncycles=6; %6;
deltaT=1*day;
nbj_buildUp=60*day;nbj_rest=20*day;nbj_inject=30*day;
nbj_idle=20*day;nbj_prod=30*day;nbj_idle1=20*day;
[schedule,TotalTime,nbuildUp,nrest,ninject,nidle,nprod,nidle1] = ...
    createCyclicScenario2( deltaT, ncycles, nbj_buildUp,...
    nbj_rest, nbj_inject,nbj_idle, nbj_prod,nbj_idle1, [W1;W2;W3;W5]);

%% Model Setup: Compositional Model with Bacterial Growth
    eosname='sw';
    model_nbs0.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_bs0.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_bs01e9.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_nbs3.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_bs3.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
    
    arg_nobact = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',false,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    arg_bact = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',true,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};          
     
    model_nbs0 = BiochemistryModel(arg_nobact{:});
    model_nbs0.outputFluxes = false;
    model_nbs0.EOSModel.msalt=0;

    model_bs0 = BiochemistryModel(arg_bact{:});
    model_bs0.outputFluxes = false;
    model_bs0.EOSModel.msalt=0;

    model_bs01e9 = BiochemistryModel(arg_bact{:});
    model_bs01e9.outputFluxes = false;
    model_bs01e9.EOSModel.msalt=0;

    model_nbs3 = BiochemistryModel(arg_nobact{:});
    model_nbs3.outputFluxes = false;
    model_nbs3.EOSModel.msalt=3;

    model_bs3 = BiochemistryModel(arg_bact{:});
    model_bs3.outputFluxes = false;
    model_bs3.EOSModel.msalt=3;

   
%% Initial Conditions
% Temperature and initial saturations
T0 = 40+273.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.7, 0.0, 0.02, 0.28];  % Initial composition: H2O, H2, CO2, C1

% Initialize state with bacterial concentration
if model_bs0.bacteriamodel
    nbact0 = 1; 
    model_bs0.nbactMax=1.e8;
    state0_bs0 = initCompositionalStateBacteria(model_bs0, P0, T0, s0, ...
        z0, nbact0,model_bs0.EOSModel);
end

if model_bs01e9.bacteriamodel
    nbact0 = 1;
    model_bs01e9.nbactMax=1.e9;
    state0_bs01e9 = initCompositionalStateBacteria(model_bs01e9, P0, T0, s0, ...
        z0, nbact0,model_bs01e9.EOSModel);
end

if ~model_nbs0.bacteriamodel
    state0_nbs0 = initCompositionalState(model_nbs0, P0, T0, s0, z0);
end
if model_bs3.bacteriamodel
    nbact0 = 1; 
    model_bs3.nbactMax=1.e8;
    state0_bs3 = initCompositionalStateBacteria(model_bs3, P0, T0, s0, ...
        z0, nbact0,model_bs3.EOSModel);
end

if ~model_nbs3.bacteriamodel
    state0_nbs3 = initCompositionalState(model_nbs3, P0, T0, s0, z0);
end


%% Run simulation
%% Pack the simulation problem with the defined components
name_nbs0='Benchmark2023AEGE_180_pack_NOBACT_6cycles';
name_bs0='Benchmark2023AEGE_180_pack_6cycles_n0_1e8';
name_bs01e9='Benchmark2023AEGE_180_pack_6cycles_n0_1e9';
name_nbs3='Benchmark2023AEGE_180_pack_NOBACT_6cycles_msalt3';
name_bs3='Benchmark2023AEGE_180_pack_6cycles_n0_1e8_msalt3';
problem_nbs0 = packSimulationProblem(state0_nbs0, model_nbs0, schedule, name_nbs0, 'NonLinearSolver', nls);
problem_bs0= packSimulationProblem(state0_bs0, model_bs0, schedule, name_bs0, 'NonLinearSolver', nls);
problem_bs01e9= packSimulationProblem(state0_bs01e9, model_bs01e9, schedule, name_bs01e9, 'NonLinearSolver', nls);
problem_nbs3 = packSimulationProblem(state0_nbs3, model_nbs3, schedule, name_nbs3, 'NonLinearSolver', nls);
problem_bs3 = packSimulationProblem(state0_bs3, model_bs3, schedule, name_bs3, 'NonLinearSolver', nls);


%% Execute the simulation of the packed problem
simulatePackedProblem(problem_nbs0); %recupere resultats si  le repertoire existe
simulatePackedProblem(problem_bs0); 
simulatePackedProblem(problem_bs01e9); 
simulatePackedProblem(problem_nbs3); 
simulatePackedProblem(problem_bs3); 

%% Get reservoir and well states
[ws_bs0,states_bs0] = getPackedSimulatorOutput(problem_bs0);
[ws_bs01e9,states_bs01e9] = getPackedSimulatorOutput(problem_bs01e9);
[ws_nbs0,states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
[ws_bs3,states_bs3] = getPackedSimulatorOutput(problem_bs3);
[ws_nbs3,states_nbs3] = getPackedSimulatorOutput(problem_nbs3);


    %% initialise unknowns %%%%%%%%%%%%%%%%%%
    namecp = model_bs0.EOSModel.getComponentNames();
    ncomp=model_bs0.EOSModel.getNumberOfComponents;

    indH2=find(strcmp(namecp,'H2'));
    indCO2= find(strcmp(namecp,'CO2'));
    indCH4= find(strcmp(namecp,'C1'));
    nT=numel(states_nbs0);

    xH2_nbs0=zeros(nT,1);
    yH2_nbs0=zeros(nT,1);
    pressure_nbs0=zeros(nT,1);
    H2_well_nbs0=zeros(nT,1);
    H2_well_bs0=zeros(nT,1);
    H2_well_bs01e9=zeros(nT,1);
    CO2_well_nbs0=zeros(nT,1);
    CO2_well_bs0=zeros(nT,1);
    CO2_well_bs01e9=zeros(nT,1);
    CH4_well_nbs0=zeros(nT,1);
    CH4_well_bs0=zeros(nT,1);
    CH4_well_bs01e9=zeros(nT,1);

    pressure_nbs3=zeros(nT,1);
    H2_well_nbs3=zeros(nT,1);
    H2_well_bs3=zeros(nT,1);
    CO2_well_nbs3=zeros(nT,1);
    CO2_well_bs3=zeros(nT,1);
    CH4_well_nbs3=zeros(nT,1);
    CH4_well_bs3=zeros(nT,1);

    xH2_bs0=zeros(nT,1);
    yH2_bs0=zeros(nT,1);
    pressure_bs0=zeros(nT,1);
    pressure_bs01e9=zeros(nT,1);
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
    totMassH2_bs01e9= zeros(nT,1);
    totMassCO2_bs01e9= zeros(nT,1);
    totMassCH4_bs01e9= zeros(nT,1);
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
        pressure_nbs0(i)=mean(states_nbs0{i}.pressure(:));
        H2_well_nbs0(i)=ws_nbs0{i}.H2.*schedule.step.val(i);
        H2_well_bs0(i)=ws_bs0{i}.H2.*schedule.step.val(i);
        H2_well_bs01e9(i)=ws_bs01e9{i}.H2.*schedule.step.val(i);
        CO2_well_nbs0(i)=ws_nbs0{i}.CO2.*schedule.step.val(i);
        CO2_well_bs0(i)=ws_bs0{i}.CO2.*schedule.step.val(i);
        CO2_well_bs01e9(i)=ws_bs01e9{i}.CO2.*schedule.step.val(i);
        CH4_well_nbs0(i)=ws_nbs0{i}.C1.*schedule.step.val(i);
        CH4_well_bs0(i)=ws_bs0{i}.C1.*schedule.step.val(i);
        CH4_well_bs01e9(i)=ws_bs01e9{i}.C1.*schedule.step.val(i);

        pressure_nbs3(i)=mean(states_nbs3{i}.pressure(:));
        H2_well_nbs3(i)=ws_nbs3{i}.H2.*schedule.step.val(i);
        H2_well_bs3(i)=ws_bs3{i}.H2.*schedule.step.val(i);
        CO2_well_nbs3(i)=ws_nbs3{i}.CO2.*schedule.step.val(i);
        CO2_well_bs3(i)=ws_bs3{i}.CO2.*schedule.step.val(i);
        CH4_well_nbs3(i)=ws_nbs3{i}.C1.*schedule.step.val(i);
        CH4_well_bs3(i)=ws_nbs3{i}.C1.*schedule.step.val(i);

        xH2_bs0(i)=max(states_bs0{i}.x(:,indH2));
        yH2_bs0(i)=max(states_bs0{i}.y(:,indH2));
        pressure_bs0(i)=mean(states_bs0{i}.pressure(:));
        pressure_bs01e9(i)=mean(states_bs01e9{i}.pressure(:));
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
        totMassH2_bs01e9(i)=sum(states_bs01e9{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs01e9(i)=sum(states_bs01e9{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_bs01e9(i)=sum(states_bs01e9{i}.FlowProps.ComponentTotalMass{indCH4});
        FractionMassH2_bs0(i)=totMassH2_bs0(i)/totMassComp_bs0(i);
        FractionMassCO2_bs0(i)=totMassCO2_bs0(i)/totMassComp_bs0(i);

        totMassH2_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indCO2});
        totMassCH4_bs3(i)=sum(states_bs3{i}.FlowProps.ComponentTotalMass{indCH4});

        FractionMassH2_bs3(i)=totMassH2_bs3(i)/totMassComp_bs3(i);
        FractionMassCO2_bs3(i)=totMassCO2_bs3(i)/totMassComp_bs3(i);
        FractionMassCH4_bs3(i)=totMassCH4_bs3(i)/totMassComp_bs3(i);  
    end

    %% Calculate H2 production efficiency  in the well
ndeb0=nbuildUp+nrest;
njcycle=ninject+nidle+nprod+nidle1;
mH2_well_injected_nbs0=sum(H2_well_nbs0(1:nbuildUp));
mH2_well_injected_nbs3=sum(H2_well_nbs3(1:nbuildUp));
mH2_well_injected_bs0=sum(H2_well_bs0(1:nbuildUp));
mH2_well_injected_bs01e9=sum(H2_well_bs01e9(1:nbuildUp));
mH2_well_injected_bs3=sum(H2_well_bs3(1:nbuildUp));
mH2_well_produced_nbs0=0.0;
mH2_well_produced_nbs3=0.0;
mH2_well_produced_bs0=0.0;
mH2_well_produced_bs01e9=0.0;
mH2_well_produced_bs3=0.0;
Efficiency_H2_well_nbs0_cycle=zeros(ncycles,1);
Efficiency_H2_well_bs0_cycle=zeros(ncycles,1);

for cycle=1:ncycles
    ndebi=ndeb0+(cycle-1)*njcycle; nj1=ndebi+ninject;
    mH2_well_injected_nbs0=mH2_well_injected_nbs0+...
        sum(H2_well_nbs0(ndebi+1:nj1));
    mH2_well_produced_nbs0=mH2_well_produced_nbs0+...
        sum(H2_well_nbs0(nj1+nidle+1:nj1+nidle+nprod));
    mH2_well_injected_nbs3=mH2_well_injected_nbs3+...
        sum(H2_well_nbs3(ndebi+1:nj1));
    mH2_well_produced_nbs3=mH2_well_produced_nbs3+...
        sum(H2_well_nbs3(nj1+nidle+1:nj1+nidle+nprod));
    mH2_well_injected_bs0=mH2_well_injected_bs0+...
        sum(H2_well_bs0(ndebi+1:nj1));
    mH2_well_produced_bs0=mH2_well_produced_bs0+...
        sum(H2_well_bs0(nj1+nidle+1:nj1+nidle+nprod));
    mH2_well_injected_bs01e9=mH2_well_injected_bs01e9+...
        sum(H2_well_bs01e9(ndebi+1:nj1));
    mH2_well_produced_bs01e9=mH2_well_produced_bs01e9+...
        sum(H2_well_bs01e9(nj1+nidle+1:nj1+nidle+nprod));
    mH2_well_injected_bs3=mH2_well_injected_bs3+...
        sum(H2_well_bs3(ndebi+1:nj1));
    mH2_well_produced_bs3=mH2_well_produced_bs3+...
        sum(H2_well_bs3(nj1+nidle+1:nj1+nidle+nprod));
    Efficiency_H2_well_nbs0_cycle(cycle)=abs(mH2_well_produced_nbs0./mH2_well_injected_nbs0).*100;
    Efficiency_H2_well_bs0_cycle(cycle)=abs(mH2_well_produced_bs0./mH2_well_injected_bs0).*100;
end
Efficiency_H2_well_nbs0=abs(mH2_well_produced_nbs0./mH2_well_injected_nbs0).*100;
Efficiency_H2_well_nbs3=abs(mH2_well_produced_nbs3./mH2_well_injected_nbs3).*100;
Efficiency_H2_well_bs0=abs(mH2_well_produced_bs0./mH2_well_injected_bs0).*100;
Efficiency_H2_well_bs01e9=abs(mH2_well_produced_bs01e9./mH2_well_injected_bs01e9).*100;
Efficiency_H2_well_bs3=abs(mH2_well_produced_bs3./mH2_well_injected_bs3).*100;
%% Display H2 production efficiency
fprintf('In the well, H2 Production Efficiency_nobact, msalt=0 : %.2f%%\n', Efficiency_H2_well_nbs0);   
fprintf('In the well, H2 Production Efficiency_nobact, msalt=3 : %.2f%%\n', Efficiency_H2_well_nbs3);   
fprintf('In the well, H2 Production Efficiency_bact, msalt=0 : %.2f%%\n', Efficiency_H2_well_bs0);   
fprintf('In the well, H2 Production Efficiency_bact, msalt=0, nbact1e8 : %.2f%%\n', Efficiency_H2_well_bs01e9);   
fprintf('In the well, H2 Production Efficiency_bact, msalt=3 : %.2f%%\n', Efficiency_H2_well_bs3);   


    %% Calculate percentage of H2 loss
    H2_loss_percentage_nbs0 = (abs(totMassH2_nbs0-totMassH2_bs0)./totMassH2_nbs0) * 100;
    H2_loss_percentage_nbs01e9 = (abs(totMassH2_nbs0-totMassH2_bs01e9)./totMassH2_nbs0) * 100;
    H2_loss_percentage_nbs3 = (abs(totMassH2_nbs3-totMassH2_bs3)./totMassH2_nbs3) * 100;
    %% Calculate percentage of CO2 loss
    CO2_loss_percentage_nbs0 = (abs(totMassCO2_nbs0-totMassCO2_bs0)./totMassCO2_nbs0) * 100;
    CO2_loss_percentage_nbs01e9 = (abs(totMassCO2_nbs0-totMassCO2_bs01e9)./totMassCO2_nbs0) * 100;
    CO2_loss_percentage_nbs3 = (abs(totMassCO2_nbs3-totMassCO2_bs3)./totMassCO2_nbs3) * 100;
    %% Calculate percentage of CH4 production
    CH4_loss_percentage_nbs0 = (abs(totMassCH4_nbs0-totMassCH4_bs0)./totMassCH4_nbs0) * 100;
    CH4_loss_percentage_nbs01e9 = (abs(totMassCH4_nbs0-totMassCH4_bs01e9)./totMassCH4_nbs0) * 100;
    CH4_loss_percentage_nbs3 = (abs(totMassCH4_nbs3-totMassCH4_bs3)./totMassCH4_nbs3) * 100;

%% Display final H2 loss
fprintf('Total H2 loss due to bacterial effects, msalt=0: %.2f%%\n',...
    H2_loss_percentage_nbs0(end));
fprintf('Total CO2 loss due to bacterial effects, msalt=0: %.2f%%\n', ...
    CO2_loss_percentage_nbs0(end));
fprintf('Total CH4 production due to bacterial effects, msalt=0: %.2f%%\n',...
    CH4_loss_percentage_nbs0(end));
fprintf('Total H2 loss due to bacterial effects, msalt=3: %.2f%%\n',...
    H2_loss_percentage_nbs3(end));
fprintf('Total CO2 loss due to bacterial effects, msalt=3: %.2f%%\n', ...
    CO2_loss_percentage_nbs3(end));
fprintf('Total CH4 production due to bacterial effects, msalt=3: %.2f%%\n',...
    CH4_loss_percentage_nbs3(end));
 fprintf('Total H2 loss due to bacterial effects, msalt=0, nbact0=1e8: %.2f%%\n',...
     H2_loss_percentage_nbs01e9(end));
fprintf('Total CO2 loss due to bacterial effects, msalt=0, nbact0=1e8: %.2f%%\n',...
    CO2_loss_percentage_nbs01e9(end));
fprintf('Total CH4 production due to bacterial effects, msalt=0, nbact0=1e8: %.2f%%\n',...
    CH4_loss_percentage_nbs01e9(end));


%% plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot pressure
f01=figure('Name','Pressure_compar_bact_msalt0_nbact0','NumberTitle','off');
f01.Position(3:4) = [900 700];
plot(1:nT,pressure_bs0./1e5,'b','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_bs01e9./1e5,'g','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_nbs0./1e5,'r--','MarkerSize',8,'LineWidth',2)
title('Mean pressure in the pure water reservoir','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'pressure (bar)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'pressure, with archae, n0=1e8','pressure, with archae,nb0=1e9','pressure, no archae'},...
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


%% plot loss H2, CO2, production C1%%%%%%%%%%%%%%%%
f20=figure('Name','H2_loss','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
%ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2 , msalt=0','H2, msalt=3'},...
    'FontSize',16,'TextColor','black','Location','west')


f21=figure('Name','H2_loss','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(1:nT,H2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_loss_percentage_nbs01e9,'g','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 16; 
legend({'H2 , msalt=0, n0=1e8','H2 , msalt=0, n0=1e9','H2, msalt=3, n0=1e9'},...
    'FontSize',16,'TextColor','black','Location','west')


f21=figure('Name','CO2_loss','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(1:nT,CO2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_percentage_nbs01e9,'g','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 16; 
legend({'CO2 , msalt=0, n0=1e8','CO2 , msalt=0, n0=1e9','CO2, msalt=3, n0=1e9'},...
    'FontSize',16,'TextColor','black','Location','west')

f21=figure('Name','C1_prod','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(1:nT,CH4_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CH4_loss_percentage_nbs01e9,'g','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CH4_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' C1 production','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'C1 production (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 16; 
legend({'C1 , msalt=0, n0=1e8','C1 , msalt=0, n0=1e9','C1, msalt=3, n0=1e9'},...
    'FontSize',16,'TextColor','black','Location','west')




f22=figure('Name','CO2_loss','NumberTitle','off');
f22.Position(3:4) = [900 700];
plot(1:nT,CO2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'CO2 , msalt=0','CO2, msalt=3'},...
    'FontSize',16,'TextColor','black','Location','west')


f23=figure('Name','CO2_loss','NumberTitle','off');
f23.Position(3:4) = [900 700];
plot(1:nT,CO2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_percentage_nbs01e9,'g','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'CO2 , msalt=0, nbact0=1e9','CO2 , msalt=0, nbact0=1e8','CO2, msalt=3, nbact0=1e9'},...
    'FontSize',16,'TextColor','black','Location','west')

f20=figure('Name','CO2_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassCO2_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassCO2_bs0,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 totmass','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'CO2 , no bact','CO2, boact'},...
    'FontSize',16,'TextColor','black','Location','west')



f20=figure('Name','C1_prod','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,CH4_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CH4_loss_percentage_nbs3,'r--','MarkerSize',7,'LineWidth',2)
title(' C1 production','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'C1 production (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'C1 , msalt=0','C1, msalt=3'},...
    'FontSize',16,'TextColor','black','Location','west')


%% plot total masses over time of: H2, CO2, production C1%%%%%%%%%%%%%%%%
f20=figure('Name','H2_total_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassH2_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassH2_bs0,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'H2 total mass, no archae','H2 total mass, archae'},...
    'FontSize',16,'TextColor','black','Location','west')


f20=figure('Name','CO2_total_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassCO2_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassCO2_bs0,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'CO2 total mass, no archae','CO2 total mass, archae'},...
    'FontSize',16,'TextColor','black','Location','west')

f20=figure('Name','C1_total_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassCH4_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassCH4_bs0,'r--','MarkerSize',7,'LineWidth',2)
title(' C1 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'C1 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'C1 total mass, no archae','C1 total mass, archae'},...
    'FontSize',16,'TextColor','black','Location','west')





%% plot well%%%%%%%%%%%%%%%%
f20=figure('Name','H2_well','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs0./schedule.step.val,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_nbs0./schedule.step.val,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2, with archae, n0=1e8','H2, no archae'},...
    'FontSize',16,'TextColor','black','Location','west')

f21=figure('Name','CO2_well','NumberTitle','off');
f21.Position(3:4) = [900 700];
plot(1:nT,CO2_well_bs0./schedule.step.val,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CO2_well_nbs0./schedule.step.val,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'CO2, with archae, n0=1e8','CO2, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')

f22=figure('Name','C1_well','NumberTitle','off');
f22.Position(3:4) = [900 700];
plot(1:nT,CH4_well_bs0./schedule.step.val,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,CH4_well_nbs0./schedule.step.val,'r--','MarkerSize',7,'LineWidth',2)
title(' C1 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'C1 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'C1, with archae, n0=1e8','C1, no archae'},...
    'FontSize',16,'TextColor','black','Location','best')


f25=figure('Name','H2_well_bact_comparsalt','NumberTitle','off');
f25.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_bs3,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate with archae','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2, msalt=0','H2, msalt=3'},...
    'FontSize',16,'TextColor','black','Location','west')


if model_bs0.bacteriamodel
    nbacteria_bs0= zeros(nT,1);
    nbacteria_bs01e9= zeros(nT,1);
    nbacteria_bs3= zeros(nT,1);
    pv=model_bs0.operators.pv;
    ncells=G.cells.num;
    for i = 1:nT
       nbacteria_bs0(i)=sum(states_bs0{i}.nbact);
       nbacteria_bs01e9(i)=sum(states_bs01e9{i}.nbact);
       nbacteria_bs3(i)=sum(states_bs3{i}.nbact);
    end

   f31=figure('Name','nbacteria','NumberTitle','off');
   f31.Position(3:4) = [900 700];
   plot(0:nT,[ncells*nbact0;nbacteria_bs0],'b','MarkerSize',7,'LineWidth',2)
   hold on
   plot(0:nT,[ncells*nbact0;nbacteria_bs01e9],'g','MarkerSize',7,'LineWidth',2)
   hold on
   plot(0:nT,[ncells*nbact0;nbacteria_bs3],'r--','MarkerSize',7,'LineWidth',2)
   title('Total methanogenic Archae population','FontSize',16,'FontWeight','bold','Color','k')
   xlabel({'time (days)'},'FontWeight','bold','Color','k')
   ylabel({'N_{archae}'},'FontWeight','bold','Color','k')
   ax = gca;
   ax.FontSize = 16;
   legend({'N_{archae}, msalt=0, n0=1e8','N_{archae}, msalt=0, nbact0=1e9',...
       'N_{archae}, msalt=3, nbact0=1e8'},'FontSize',16,'TextColor','black',...
    'Location','best')




f33=figure('Name','nbacteria_t110','NumberTitle','off');
f33.Position(3:4) = [900 700];
plotCellData(G,states_bs0{110}.nbact);
colorbar; 
axis equal
axis ([0 Lx  0 Ly depth_res depth_res+Lz])
view(0,-90)
title('Methanogenic Archae population, 110 days','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'x (m)'},'FontWeight','bold','Color','k')
ylabel({'y (m)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 

f33=figure('Name','xH2_t110','NumberTitle','off');
f33.Position(3:4) = [900 700];
plotCellData(G,states_bs0{110}.x(:,indH2));
colorbar; 
axis equal
axis ([0 Lx  0 Ly depth_res depth_res+Lz])
view(0,-90)
title('H2 molar fraction in the liquid phase, 110 days','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'x (m)'},'FontWeight','bold','Color','k')
ylabel({'y (m)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 

end