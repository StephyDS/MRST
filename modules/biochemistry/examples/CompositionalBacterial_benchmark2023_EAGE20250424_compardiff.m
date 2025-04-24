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
nobact=true;%false;%

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

ncycles=2; %6;
deltaT=1*day;
nbj_buildUp=60*day;nbj_rest=20*day;nbj_inject=30*day;
nbj_idle=20*day;nbj_prod=30*day;nbj_idle1=20*day;
[schedule,TotalTime,nbuildUp,nrest,ninject,nidle,nprod,nidle1] = ...
    createCyclicScenario2( deltaT, ncycles, nbj_buildUp,...
    nbj_rest, nbj_inject,nbj_idle, nbj_prod,nbj_idle1, [W1;W2;W3;W5]);

%% Model Setup: Compositional Model with Bacterial Growth
    eosname='sw';
    model_nbs0.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_nbs0_moldiffmq.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    model_bs0.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
    
    arg_nobact = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',false,...
        'bDiffusionEffect', false,'moleculardiffusion',false,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    arg_nobact_moldiff = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',false,...
        'bDiffusionEffect', false,'moleculardiffusion',true,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};
    arg_bact = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',true,...
        'bDiffusionEffect', false,'moleculardiffusion',true,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};          
     
    model_nbs0 = BiochemistryModel(arg_nobact{:});
    model_nbs0.outputFluxes = false;
    model_nbs0.EOSModel.msalt=0;

    model_nbs0_moldiffmq = BiochemistryModel(arg_nobact{:});
    model_nbs0_moldiffmq.outputFluxes = false;
    model_nbs0_moldiffmq.EOSModel.msalt=0;

    model_bs0 = BiochemistryModel(arg_bact{:});
    model_bs0.outputFluxes = false;
    model_bs0.EOSModel.msalt=0;

     
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

if ~model_nbs0.bacteriamodel
    state0_nbs0 = initCompositionalState(model_nbs0, P0, T0, s0, z0);
    state0_nbs0_moldiffmq = initCompositionalState(model_nbs0_moldiffmq, P0, T0, s0, z0);
end



%% Run simulation
%% Pack the simulation problem with the defined components
name_nbs0='Benchmark2023AEGE_180_pack_NOBACT_2cycles';
name_nbs0_moldiffmq='Benchmark2023AEGE_180_pack_NOBACT_2cycles_MolDiff_MQ';
name_bs0='Benchmark2023AEGE_180_pack_2cycles_n0_1e8_MolDiff_MQ';
problem_nbs0 = packSimulationProblem(state0_nbs0, model_nbs0, schedule, name_nbs0, 'NonLinearSolver', nls);
problem_nbs0_moldiffmq = packSimulationProblem(state0_nbs0_moldiffmq,...
    model_nbs0_moldiffmq, schedule, name_nbs0_moldiffmq, 'NonLinearSolver', nls);
problem_bs0= packSimulationProblem(state0_bs0, model_bs0, schedule, name_bs0, 'NonLinearSolver', nls);

%% Execute the simulation of the packed problem
simulatePackedProblem(problem_nbs0); %recupere resultats si  le repertoire existe
simulatePackedProblem(problem_nbs0_moldiffmq); %recupere resultats si  le repertoire existe
simulatePackedProblem(problem_bs0); 

%% Get reservoir and well states
[ws_bs0,states_bs0] = getPackedSimulatorOutput(problem_bs0);
[ws_nbs0,states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
[ws_nbs0_moldiffmq,states_nbs0_moldiffmq] = getPackedSimulatorOutput(problem_nbs0_moldiffmq);


    %% initialise unknowns %%%%%%%%%%%%%%%%%%
    namecp = model_nbs0.EOSModel.getComponentNames();
    ncomp=model_nbs0.EOSModel.getNumberOfComponents;

    indH2=find(strcmp(namecp,'H2'));
    indCO2= find(strcmp(namecp,'CO2'));
    indCH4= find(strcmp(namecp,'C1'));
    nT=numel(states_nbs0);

    xH2_nbs0=zeros(nT,1);
    yH2_nbs0=zeros(nT,1);
    pressure_nbs0=zeros(nT,1);
    xH2_nbs0_moldiffmq=zeros(nT,1);
    yH2_nbs0_moldiffmq=zeros(nT,1);
    pressure_nbs0_moldiffmq=zeros(nT,1);
    xH2_bs0=zeros(nT,1);
    yH2_bs0=zeros(nT,1);
    pressure_bs0=zeros(nT,1);
   
    H2_well_nbs0=zeros(nT,1);
    H2_well_bs0=zeros(nT,1);
    H2_well_nbs0_moldiffmq=zeros(nT,1);
   
   
    
    totMassH2_nbs0= zeros(nT,1);
    totMassCO2_nbs0= zeros(nT,1);
    totMassH2_bs0= zeros(nT,1);
    totMassCO2_bs0= zeros(nT,1);
   
    totMassH2_nbs0_moldiffmq= zeros(nT,1);
    totMassCO2_nbs0_moldiffmq= zeros(nT,1);
   
       for i = 1:nT
        xH2_nbs0(i)=max(states_nbs0{i}.x(:,indH2));
        yH2_nbs0(i)=max(states_nbs0{i}.y(:,indH2));
        xH2_nbs0_moldiffmq(i)=max(states_nbs0_moldiffmq{i}.x(:,indH2));
        yH2_nbs0_moldiffmq(i)=max(states_nbs0_moldiffmq{i}.y(:,indH2));
        pressure_nbs0_moldiffmq(i)=mean(states_nbs0_moldiffmq{i}.pressure(:));
        pressure_nbs0(i)=mean(states_nbs0{i}.pressure(:));
        pressure_bs0(i)=mean(states_bs0{i}.pressure(:));        

        H2_well_nbs0(i)=ws_nbs0{i}.H2.*schedule.step.val(i);
        H2_well_bs0(i)=ws_bs0{i}.H2.*schedule.step.val(i);
        H2_well_nbs0_moldiffmq(i)=ws_nbs0_moldiffmq{i}.H2.*schedule.step.val(i);
        
        
      end

    %% mass calculations
    for i = 1:nT
        totMassH2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indCO2});
               
        totMassH2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indCO2});
        
        totMassH2_nbs0_moldiffmq(i)=sum(states_nbs0_moldiffmq{i}.FlowProps.ComponentTotalMass{indH2});
        totMassCO2_nbs0_moldiffmq(i)=sum(states_nbs0_moldiffmq{i}.FlowProps.ComponentTotalMass{indCO2});
     end

    %% Calculate H2 production efficiency  in the well
ndeb0=nbuildUp+nrest;
njcycle=ninject+nidle+nprod+nidle1;
mH2_well_injected_nbs0=sum(H2_well_nbs0(1:nbuildUp));
mH2_well_injected_bs0=sum(H2_well_bs0(1:nbuildUp));
mH2_well_injected_nbs0_moldiffmq=sum(H2_well_nbs0_moldiffmq(1:nbuildUp));
mH2_well_produced_nbs0=0.0;
mH2_well_produced_bs0=0.0;
mH2_well_produced_nbs0_moldiffmq=0.0;
Efficiency_H2_well_nbs0_cycle=zeros(ncycles,1);
Efficiency_H2_well_bs0_cycle=zeros(ncycles,1);
Efficiency_H2_well_nbs0_moldiffmq_cycle=zeros(ncycles,1);

for cycle=1:ncycles
    ndebi=ndeb0+(cycle-1)*njcycle; nj1=ndebi+ninject;
    mH2_well_injected_nbs0=mH2_well_injected_nbs0+...
        sum(H2_well_nbs0(ndebi+1:nj1));
    mH2_well_produced_nbs0=mH2_well_produced_nbs0+...
        sum(H2_well_nbs0(nj1+nidle+1:nj1+nidle+nprod));
   
    mH2_well_injected_bs0=mH2_well_injected_bs0+...
        sum(H2_well_bs0(ndebi+1:nj1));
    mH2_well_produced_bs0=mH2_well_produced_bs0+...
        sum(H2_well_bs0(nj1+nidle+1:nj1+nidle+nprod));
   
    mH2_well_injected_nbs0_moldiffmq=mH2_well_injected_nbs0_moldiffmq+...
        sum(H2_well_nbs0_moldiffmq(ndebi+1:nj1));
    mH2_well_produced_nbs0_moldiffmq=mH2_well_produced_nbs0_moldiffmq+...
        sum(H2_well_nbs0_moldiffmq(nj1+nidle+1:nj1+nidle+nprod));

    Efficiency_H2_well_nbs0_cycle(cycle)=abs(mH2_well_produced_nbs0./mH2_well_injected_nbs0).*100;
    Efficiency_H2_well_bs0_cycle(cycle)=abs(mH2_well_produced_bs0./mH2_well_injected_bs0).*100; 
    Efficiency_H2_well_nbs0_moldiffmq_cycle(cycle)=...
        abs(mH2_well_produced_nbs0_moldiffmq./mH2_well_injected_nbs0_moldiffmq).*100;


end
%% Display H2 production efficiency
fprintf('In the well, H2 Production Efficiency_nobact, no molecular diffusion : %.2f%%\n', Efficiency_H2_well_nbs0_cycle(end));   
fprintf('In the well, H2 Production Efficiency_nobact, molecular diffusion MQ: %.2f%%\n', Efficiency_H2_well_nbs0_moldiffmq_cycle(end));   
fprintf('In the well, H2 Production Efficiency_bact, molecular diffusion MQ : %.2f%%\n', Efficiency_H2_well_bs0_cycle(end));   


    %% Calculate percentage of H2 loss
    H2_loss_percentage_nbs0 = (abs(totMassH2_nbs0-totMassH2_bs0)./totMassH2_nbs0) * 100;
    H2_loss_percentage_nbs0MQ = (abs(totMassH2_nbs0-totMassH2_nbs0_moldiffmq)./totMassH2_nbs0) * 100;
   %% Calculate percentage of CO2 loss
    CO2_loss_percentage_nbs0 = (abs(totMassCO2_nbs0-totMassCO2_bs0)./totMassCO2_nbs0) * 100;
    CO2_loss_percentage_nbs0MQ = (abs(totMassCO2_nbs0-totMassCO2_nbs0_moldiffmq)./totMassCO2_nbs0) * 100;
     
%% Display final H2 loss
fprintf('Total H2 loss due to molecular diffusion AND bacteria: %.2f%%\n',...
    H2_loss_percentage_nbs0(end));
fprintf('Total H2 loss due to molecular diffusion: %.2f%%\n',...
    H2_loss_percentage_nbs0MQ(end));
fprintf('Total CO2 loss due to molecular diffusion AND bacteria: %.2f%%\n', ...
    CO2_loss_percentage_nbs0(end));
fprintf('Total CO2 loss due to molecular diffusion: %.2f%%\n', ...
    CO2_loss_percentage_nbs0MQ(end));


%% plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot pressure
f01=figure('Name','Pressure_compar_diff','NumberTitle','off');
f01.Position(3:4) = [900 700];
plot(1:nT,pressure_nbs0./1e5,'b','MarkerSize',7,'LineWidth',2)
hold on
plot(1:nT,pressure_bs0./1e5,'r--','MarkerSize',8,'LineWidth',2)
hold on
plot(1:nT,pressure_nbs0_moldiffmq./1e5,'g--','MarkerSize',8,'LineWidth',2)
title('Mean pressure in the pure water reservoir','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'pressure (bar)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'pressure, no molecular diffusion, no bact','pressure, molecular diffusion and bacteria','pressure, molecular diffusion, no bact'},...
    'FontSize',16,'TextColor','black',...
    'Location','best')

%% plot efficiency
% f01=figure('Name','efficiency_msalt0','NumberTitle','off');
% f01.Position(3:4) = [900 700];
% plot(1:6,Efficiency_H2_well_nbs0_cycle(1:6),'r--','MarkerSize',7,'LineWidth',2)
% hold on
% plot(1:6,Efficiency_H2_well_nbs0_moldiff_cycle(1:6),'b','MarkerSize',7,'LineWidth',2)
% title('Efficiency production of H2','FontSize',16,'FontWeight','bold','Color','k')
% xlabel({'cycles'},'FontWeight','bold','Color','k')
% ylabel({'efficiency H2(%)'},'FontWeight','bold','Color','k')
% ax = gca;
% ax.FontSize = 16; 
% legend({'efficiency, no archae','efficiency, with archae,nb0=1e8'},...
%     'FontSize',16,'TextColor','black',...
%     'Location','best')

%% plot loss H2, CO2, production C1%%%%%%%%%%%%%%%%
f20=figure('Name','H2_loss','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_loss_percentage_nbs0MQ,'r--','MarkerSize',7,'LineWidth',2)
title(' H2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.FontSize = 16; 
legend({'H2 , constant molecular diffusion and bacteria','H2 molecular diffusion no bact'},...
    'FontSize',16,'TextColor','black','Location','west')



%% plot total masses over time of: H2, CO2, production C1%%%%%%%%%%%%%%%%
f20=figure('Name','H2_total_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassH2_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassH2_bs0,'r--','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassH2_nbs0_moldiffmq,'g--','MarkerSize',7,'LineWidth',2)
title(' H2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'H2 total mass, no molecular diffusion, no bact',...
    'H2 total mass, molecular diffusion and bact','H2 total mass, molecular diffusion, no bact'},...
    'FontSize',16,'TextColor','black','Location','west')


f20=figure('Name','CO2_total_mass','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,totMassCO2_nbs0,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,totMassCO2_nbs0_moldiff,'r--','MarkerSize',7,'LineWidth',2)
title(' CO2 total mass over time, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'CO2 mass (kg)'},'FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 16; 
legend({'CO2 total mass, molecular diffusion','CO2 total mass, molecular diffusion'},...
    'FontSize',16,'TextColor','black','Location','west')


%% plot well%%%%%%%%%%%%%%%%
f20=figure('Name','H2_well','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_well_nbs0./schedule.step.val,'b','MarkerSize',7,'LineWidth',2)
hold on;
plot(1:nT,H2_well_nbs0_moldiffmq./schedule.step.val,'g--','MarkerSize',7,'LineWidth',2)
title(' H2 well rate, no salt','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2, molecular diffusion and bacteria','H2, molecular diffusion, no bacteria'},...
    'FontSize',16,'TextColor','black','Location','west')

