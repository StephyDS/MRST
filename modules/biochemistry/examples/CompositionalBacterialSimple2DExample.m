clear; clc;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%====Geometrie, maillage et proprietes roche
[nx,ny,nz] = deal(20,1,20);
[Lx,H]=deal(200,100);
dims = [nx, ny,nz];
pdims = [Lx, 1, H];
dz=H/nz;
G = cartGrid(dims, pdims);
depth_res=1100; %profondeur haut reservoir
G.nodes.coords(:,3) = G.nodes.coords(:,3) + depth_res;

G = computeGeometry(G);
rock = makeRock(G, 30*milli*darcy, 0.2);

[i1,i2,i3,i4]=deal(floor(0.05*nz),floor(0.4*nz),floor(0.5*nz),floor(0.7*nz));
for i=1:i1
    I = (1:nx).'; 
    J = ones(nx,1); %repmat(1, [nx, 1]);
    K = repmat(i, [nx, 1]);
    ind = sub2ind(G.cartDims, I, J, K);
    rock.perm(ind) = 1e-19;
    rock.poro(ind) = 0.001;
end
for i=i2+1:i3
    I = (1:nx).'; 
    J = ones(nx,1);%repmat(1, [nx, 1]);
    K = repmat(i, [nx, 1]);
    ind = sub2ind(G.cartDims, I, J, K);
    rock.perm(ind) = 1e-15;
    rock.poro(ind) = 0.05;
end
 figure;
 plotCellData(G,rock.poro,'EdgeAlpha',0.2)
 colorbar
 view(3)




%%=========initialisation proprietes du fluide
%===Compositional fluid model (initialization with the CoolProp library)=
compFluid =TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide',...
    'Nitrogen','Methane'}, {'H2O','H2','CO2','N2','CH4'});

%Brine-Gas (H2)
[rhow,rhog]=deal(1000* kilogram/meter^3,8.1688* kilogram/meter^3); %density kilogram/meter^3;
[viscow,viscog]=deal(1.0*centi*poise,0.0094234*centi*poise);%viscosity
[cfw,cfg]=deal(0,8.1533e-3/barsa); %compressibility
 
[srw,src]=deal(0.0,0.0);
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
[Pmaxz,Pref1,Pminz,Pe]=deal(95*barsa,114*barsa,120*barsa,0.1*barsa); %pressions

%initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
fluid=initSimpleADIFluid('phases', 'OG', 'mu',[viscow,viscog],...
                         'rho',[rhow,rhog],'pRef',Pref1,...
                         'c',[cfw,cfg],'n',[2,2],'smin',[srw,src]);

% Pression capillaire
pcOG = @(so) Pe * so.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1-sg-srw)./(1-srw), 1e-5)); %@@


T = 100*day;
pv=sum(poreVolume(G,rock))/T;
rate = 100*pv;
niter = 30;

%%==============Conditions aux limites et puits============
bc = [];
%puit d'injection
W = [];
W = verticalWell(W, G, rock,1,1,nz, 'comp_i', [0, 1],'Radius',0.5,...
    'name', 'Injector', 'type', 'rate','Val',rate, 'sign', 1);
W(1).components = [0.0 0.95 0.05 0.0 0.0];
% for i = 1:numel(W)
%     W(i).components = info.injection;
% end
%%==============model compositionnal================
arg = {G, rock, fluid, compFluid,...
    'water', false, 'oil', true, 'gas', true,... % water-oil system
	'bacteriamodel', true,'diffusioneffect',false,'liquidPhase', 'O',...
    'vaporPhase', 'G'}; % water=liquid, gas=vapor
model = BioChemsitryModel(arg{:});
model.outputFluxes = false;
%===Conditions initiales=====================
T0=317.5;
s0= [0.8 0.2]; %initial saturations  Sw=1
z0 = [0.8,0.0,0.006,0.018,0.176]; %initial composition: H2O,H2,CO2,N2,CH4.

%===Bacteria model===========================
if model.bacteriamodel
    nbact0=10^6;
    state0 = initCompositionalStateBacteria(model,Phydro0,T0,s0,z0,nbact0);
else
    state0 = initCompositionalState(model, Phydro0, T0, s0, z0);
end


%===Ajout d'un terme source====================
src=[];

%===Resolution pression/transport========================
deltaT = T/niter;
schedule = simpleSchedule(repmat(deltaT,1,niter),'bc', bc,'src', src,'W',W);
nls = NonLinearSolver('useRelaxation', true);
[~,states,report] = simulateScheduleAD(state0, model, schedule);

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