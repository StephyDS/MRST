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
MolecDiffus=false;%true;%
Dispers=true;%false;%
BactDiffus=false;%true;%
BactChemotaxis=false;%true;%
%% ============Grid and Rock Properties=====================
% Define grid dimensions and physical dimensions
%[nx, ny, nz] = deal(61,61,10);  % Grid cells in x, y, z directions
%[nx, ny, nz] = deal(31,31,8);  % Grid cells in x, y, z directions
[nx, ny, nz] = deal(11,11,8);  % Grid cells in x, y, z directions
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


%% Wells and Boundary Conditions
% Initialize wells
rate = 1e6*meter^3/day; 

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
%nls = NonLinearSolver('useRelaxation', true);
nls = NonLinearSolver('useRelaxation', true,'verbose',true);
%nls = NonLinearSolver();

ncycles=6; %6;
deltaT=1*day;
nbj_buildUp=60*day;nbj_rest=20*day;nbj_inject=30*day;
nbj_idle=20*day;nbj_prod=30*day;nbj_idle1=20*day;
[schedule,TotalTime,nbuildUp,nrest,ninject,nidle,nprod,nidle1] = ...
    createCyclicScenario2( deltaT, ncycles, nbj_buildUp,...
    nbj_rest, nbj_inject,nbj_idle, nbj_prod,nbj_idle1, [W1;W2;W3;W5]);

%% Model Setup: Compositional Model with Bacterial Growth
    eosname='sw';
    model.EOSModel =SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
    diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
    mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
    if nobact
        arg = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',false,...
        'bDiffusionEffect', false,'moleculardiffusion',MolecDiffus,...
        'dispersion',Dispers,'liquidPhase', 'O', 'vaporPhase', 'G'};
    else
        arg = {G, rock, fluid, compFluid,true,diagonal_backend,...
        'water', false, 'oil', true, 'gas', true,'bacteriamodel',true,...
        'bDiffusionEffect', BactDiffus,'moleculardiffusion',MolecDiffus,...
        'dispersion',Dispers,'chemotaxisEffect',BactChemotaxis,...
        'liquidPhase', 'O', 'vaporPhase', 'G'};          
    end
    model = BiochemistryModel(arg{:});
    model.outputFluxes = false;
    model.EOSModel.msalt=0;

    lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
    nls.LinearSolver = lsolve;                                     % Assign the linear solver to the nonlinear solver

   
%% Initial Conditions
% Temperature and initial saturations
T0 = 40+273.15;                % Initial temperature (K)
s0 = [0.2, 0.8];           % Initial saturations (Sw,Sg)
z0 = [0.7, 0.0, 0.02, 0.28];  % Initial composition: H2O, H2, CO2, C1

% Initialize state with bacterial concentration
if model.bacteriamodel
    nbact0 = 1; 
    model.nbactMax=5.e8;
    state0 = initCompositionalStateBacteria(model, P0, T0, s0, ...
        z0, nbact0,model.EOSModel);
else
    state0 = initCompositionalState(model, P0, T0, s0, z0);
end


%% Run simulation
%% Pack the simulation problem with the defined components
%name_nbs0='Benchmark2023AEGE_180_pack_NOBACT_6cycles_msalt3';
name_nbs0='Benchmark2023AEGE_180_pack_NOBACT_6cycles_gmres_DISP_COARSE';
problem_nbs0 = packSimulationProblem(state0, model, schedule, name_nbs0, 'NonLinearSolver', nls);
   
if nobact
    %% Execute the simulation of the packed problem
    simulatePackedProblem(problem_nbs0, 'restartStep',1);
    %% Get reservoir and well states
    [ws_nbs0,states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
else
    %name='Benchmark2023AEGE_180_pack_6cycles_n0_1e8_msalt3';
    name='Benchmark2023AEGE_180_pack_6cycles_n0_1e8_BactDiff_gmres_COARSE';
    problem_bs0 = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls); 
    %% Execute the simulation of the packed problem
      %[wellSols,states,report]= simulateScheduleAD(state0, model, schedule,...
       % 'nonlinearsolver',nls);
    simulatePackedProblem(problem_bs0, 'restartStep',1); %execute la simulation depuis le pas 1
    simulatePackedProblem(problem_nbs0); %recupere resultats si  le repertoire existe
    %% Get reservoir and well states
    [ws_bs0,states_bs0] = getPackedSimulatorOutput(problem_bs0);
    [ws_nbs0,states_nbs0] = getPackedSimulatorOutput(problem_nbs0);
end




% %% Plotting Results and estimate bacteria effects
% %% initialise unknowns %%%%%%%%%%%%%%%%%%
namecp = model.EOSModel.getComponentNames();
ncomp=model.EOSModel.getNumberOfComponents;
nT=numel(states_nbs0);
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
indC1= find(strcmp(namecp,'C1'));


xH2_nbs0=zeros(nT,1);
yH2_nbs0=zeros(nT,1);
xCO2_nbs0= zeros(nT,1);
yCO2_nbs0= zeros(nT,1);
yC1_nbs0= zeros(nT,1);
pressure_nbs0=zeros(nT,1);
H2_well_nbs0=zeros(nT,1);
CO2_well_nbs0=zeros(nT,1);
C1_well_nbs0=zeros(nT,1);
totMassH2_nbs0= zeros(nT,1);
totMassCO2_nbs0= zeros(nT,1);
totMassC1_nbs0= zeros(nT,1);
FractionMassCO2_nbs0= zeros(nT,1);
FractionMassH2_nbs0= zeros(nT,1);
FractionMassC1_nbs0= zeros(nT,1);   
totMassComp_nbs0= zeros(nT,1);

if ~nobact
    xH2_bs0=zeros(nT,1);
    yH2_bs0=zeros(nT,1);
    xCO2_bs0= zeros(nT,1);
    yCO2_bs0= zeros(nT,1);
    yC1_bs0= zeros(nT,1);
    pressure_bs0=zeros(nT,1);
    H2_well_bs0=zeros(nT,1);
    CO2_well_bs0=zeros(nT,1);
    C1_well_bs0=zeros(nT,1);
    totMassH2_bs0= zeros(nT,1);
    totMassCO2_bs0= zeros(nT,1);
    totMassC1_bs0= zeros(nT,1);
    FractionMassCO2_bs0= zeros(nT,1);
    FractionMassH2_bs0= zeros(nT,1);
    FractionMassC1_bs0= zeros(nT,1);
    totMassComp_bs0= zeros(nT,1);
end

for i = 1:nT
    xH2_nbs0(i)=max(states_nbs0{i}.x(:,indH2));
    yH2_nbs0(i)=max(states_nbs0{i}.y(:,indH2));
    yC1_nbs0(i)=max(states_nbs0{i}.y(:,indC1));
    xCO2_nbs0(i)=max(states_nbs0{i}.x(:,indCO2));
    yCO2_nbs0(i)=max(states_nbs0{i}.y(:,indCO2));
    pressure_nbs0(i)=mean(states_nbs0{i}.pressure(:));
    H2_well_nbs0(i)=ws_nbs0{i}.H2;%kg/s
    CO2_well_nbs0(i)=ws_nbs0{i}.CO2;
    C1_well_nbs0(i)=ws_nbs0{i}.C1;
end

if ~nobact
    for i = 1:nT
        xH2_bs0(i)=max(states_bs0{i}.x(:,indH2));
        yH2_bs0(i)=max(states_bs0{i}.y(:,indH2));
        yC1_bs0(i)=max(states_bs0{i}.y(:,indC1));
        xCO2_bs0(i)=max(states_bs0{i}.x(:,indCO2));
        yCO2_bs0(i)=max(states_bs0{i}.y(:,indCO2));
        pressure_bs0(i)=mean(states_bs0{i}.pressure(:));
        H2_well_bs0(i)=ws_bs0{i}.H2;
        CO2_well_bs0(i)=ws_bs0{i}.CO2;
        C1_well_bs0(i)=ws_bs0{i}.C1;
    end
end

%% mass calculations
for i = 1:nT
    for j=1:ncomp
        totMassComp_nbs0(i)=totMassComp_nbs0(i)+sum(states_nbs0{i}.FlowProps.ComponentTotalMass{j});
    end
    totMassH2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indH2});
    totMassCO2_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indCO2});
    totMassC1_nbs0(i)=sum(states_nbs0{i}.FlowProps.ComponentTotalMass{indC1});
    FractionMassH2_nbs0(i)=totMassH2_nbs0(i)/totMassComp_nbs0(i);
    FractionMassCO2_nbs0(i)=totMassCO2_nbs0(i)/totMassComp_nbs0(i);
    FractionMassC1_nbs0(i)=totMassC1_nbs0(i)/totMassComp_nbs0(i);
end

if ~nobact
    for i=1:nT
       for j=1:ncomp
           totMassComp_bs0(i)=totMassComp_bs0(i)+sum(states_bs0{i}.FlowProps.ComponentTotalMass{j});
       end
       totMassH2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indH2});
       totMassCO2_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indCO2});
       totMassC1_bs0(i)=sum(states_bs0{i}.FlowProps.ComponentTotalMass{indC1});
       FractionMassH2_bs0(i)=totMassH2_bs0(i)/totMassComp_bs0(i);
       FractionMassCO2_bs0(i)=totMassCO2_bs0(i)/totMassComp_bs0(i);
       FractionMassC1_bs0(i)=totMassC1_bs0(i)/totMassComp_bs0(i);
    end
end

% %% Calculate H2 production efficiency in the well
ndeb0=nbuildUp+nrest;
njcycle=ninject+nidle+nprod+nidle1;
mH2_well_injected_nbs0=sum(H2_well_nbs0(1:nbuildUp));
mH2_well_produced_nbs0=0.0;
for cycle=1:ncycles
    ndebi=ndeb0+(cycle-1)*njcycle; nj1=ndebi+ninject;
    mH2_well_injected_nbs0=mH2_well_injected_nbs0+...
        sum(H2_well_nbs0(ndebi+1:nj1));
    mH2_well_produced_nbs0=mH2_well_produced_nbs0+...
        sum(H2_well_nbs0(nj1+nidle+1:nj1+nidle+nprod));
end
Efficiency_H2_well_nbs0=abs(mH2_well_produced_nbs0./mH2_well_injected_nbs0).*100;
%% Display H2 production efficiency
fprintf('In the well, H2 Production Efficiency_nobact, msalt=0 : %.2f%%\n', Efficiency_H2_well_nbs0);   

if ~nobact
    %% Calculate H2 production efficiency in the well 
    ndeb0=nbuildUp+nrest;
    njcycle=ninject+nidle+nprod+nidle1;
    mH2_well_injected_bs0=sum(H2_well_bs0(1:nbuildUp));
    mH2_well_produced_bs0=0.0;
    for cycle=1:ncycles
        ndebi=ndeb0+(cycle-1)*njcycle; nj1=ndebi+ninject;
        mH2_well_injected_bs0=mH2_well_injected_bs0+...
           sum(H2_well_bs0(ndebi+1:nj1));
        mH2_well_produced_bs0=mH2_well_produced_bs0+...
           sum(H2_well_bs0(nj1+nidle+1:nj1+nidle+nprod));
    end
    Efficiency_H2_well_bs0=abs(mH2_well_produced_bs0./mH2_well_injected_bs0).*100;
    %% Display H2 production efficiency
    fprintf('In the well, H2 Production Efficiency_bact, msalt=0 : %.2f%%\n', Efficiency_H2_well_bs0);
    %% Calculate percentage of H2 loss
    H2_loss_percentage_nbs0 = (abs(totMassH2_nbs0-totMassH2_bs0)./totMassH2_nbs0) * 100;
    %% Calculate percentage of CO2 loss
    CO2_loss_percentage_nbs0 = (abs(totMassCO2_nbs0-totMassCO2_bs0)./totMassCO2_nbs0) * 100;
    %% Calculate percentage of C1 production
    C1_loss_percentage_nbs0 = (abs(totMassC1_nbs0-totMassC1_bs0)./totMassC1_nbs0) * 100;
    %% Display final H2 loss
    fprintf('Total H2 loss due to bacterial effects, msalt=0: %.2f%%\n', H2_loss_percentage_nbs0(end));
    fprintf('Total CO2 loss due to bacterial effects, msalt=0: %.2f%%\n', CO2_loss_percentage_nbs0(end));
    fprintf('Total C1 production due to bacterial effects, msalt=0: %.2f%%\n', C1_loss_percentage_nbs0(end));



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

f20=figure('Name','H2_loss','NumberTitle','off');
f20.Position(3:4) = [900 700];
plot(1:nT,H2_loss_percentage_nbs0,'b','MarkerSize',7,'LineWidth',2)
title(' H2 loss','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 loss (%)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
%ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2 , msalt=0'},...
    'FontSize',16,'TextColor','black','Location','west')


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



f25=figure('Name','H2_well_bact','NumberTitle','off');
f25.Position(3:4) = [900 700];
plot(1:nT,H2_well_bs0,'b','MarkerSize',7,'LineWidth',2)
title(' H2 well rate with archae','FontSize',16,'FontWeight','bold','Color','k')
xlabel({'time (days)'},'FontWeight','bold','Color','k')
ylabel({'H2 Production (kg/day)'},'FontWeight','bold','Color','k')
ax = gca;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.FontSize = 16; 
legend({'H2, msalt=0'},...
    'FontSize',16,'TextColor','black','Location','west')







 nbacteria_bs0= zeros(nT,1);
    pv=model.operators.pv;
    ncells=G.cells.num;
    for i = 1:nT
       nbacteria_bs0(i)=sum(states_bs0{i}.nbact);
    end

   f31=figure('Name','nbacteria','NumberTitle','off');
   f31.Position(3:4) = [900 700];
   plot(0:nT,[ncells*nbact0;nbacteria_bs0],'b','MarkerSize',7,'LineWidth',2)
   title('Total methanogenic Archae population','FontSize',16,'FontWeight','bold','Color','k')
   xlabel({'time (days)'},'FontWeight','bold','Color','k')
   ylabel({'N_{archae}'},'FontWeight','bold','Color','k')
   ax = gca;
   ax.FontSize = 16; 
   legend({'N_{archae}, msalt=0'},'FontSize',16,'TextColor','black',...
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


