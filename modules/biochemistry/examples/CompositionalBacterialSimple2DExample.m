%% MRST Simulation for Hydrogen Storage with Bacterial Growth Model
% Description: This script uses MRST to model gas injection into a porous medium,
% incorporating compositional fluid properties, bacterial mono modal.

% Clear workspace and initialize MRST modules
clear; clc;
mrstModule add compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on

%% Grid and Rock Properties
% Define grid dimensions and physical dimensions
[nx, ny, nz] = deal(20, 1, 20);  % Grid cells in x, y, z directions
[Lx, H] = deal(200, 100);         % Physical dimensions in meters
dims = [nx, ny, nz];
pdims = [Lx, 1, H];
dz = H / nz;                      % Cell height

% Create grid and shift vertically by reservoir depth
G = cartGrid(dims, pdims);
depth_res = 1100;                 % Reservoir depth in meters
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + depth_res;
G = computeGeometry(G);

% Define rock properties
rock = makeRock(G, 30 * milli * darcy, 0.2);  % Default permeability and porosity

% Custom permeability and porosity layers
[i1, i2, i3, i4] = deal(floor(0.05 * nz), floor(0.4 * nz), floor(0.5 * nz), floor(0.7 * nz));
for i = 1:i1
    % Define low permeability zone near the top of the reservoir
    ind = sub2ind(G.cartDims, (1:nx).', ones(nx, 1), repmat(i, [nx, 1]));
    rock.perm(ind) = 1e-19;
    rock.poro(ind) = 0.001;
end
for i = i2+1:i3
    % Define mid-permeability zone in the middle of the reservoir
    ind = sub2ind(G.cartDims, (1:nx).', ones(nx, 1), repmat(i, [nx, 1]));
    rock.perm(ind) = 1e-15;
    rock.poro(ind) = 0.05;
end

% Plot porosity distribution
figure;
plotCellData(G, rock.poro, 'EdgeAlpha', 0.2);
colorbar;
view(3);
title('Porosity Distribution');

%% Fluid Properties Initialization
% Define compositional fluid model (with CoolProp library support)
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Nitrogen', 'Methane'}, ...
                                      {'H2O', 'H2', 'CO2', 'N2', 'CH4'});

% Fluid density and viscosity (kg/m^3 and cP)
[rhow, rhog] = deal(1000 * kilogram / meter^3, 8.1688 * kilogram / meter^3);
[viscow, viscog] = deal(1.0 * centi * poise, 0.0094234 * centi * poise);

% Compressibility (per bar)
[cfw, cfg] = deal(0, 8.1533e-3 / barsa);

% Relative permeability and initial saturations
[srw, src] = deal(0.0, 0.0);
fluid = initSimpleADIFluid('phases', 'OG', 'mu', [viscow, viscog], ...
                           'rho', [rhow, rhog], 'pRef', 114 * barsa, ...
                           'c', [cfw, cfg], 'n', [2, 2], 'smin', [srw, src]);

% Capillary pressure function
Pe = 0.1 * barsa;
pcOG = @(so) Pe * so.^(-1/2);
fluid.pcOG = @(sg) pcOG(max((1 - sg - srw) / (1 - srw), 1e-5));

%% Simulation Parameters
% Set total time, pore volume, and injection rate
T = 100 * day;
pv = sum(poreVolume(G, rock)) / T;
rate = 100 * pv;
niter = 30;

%% Wells and Boundary Conditions
% Initialize wells
W = [];
% Injection well parameters
W = verticalWell(W, G, rock, 1, 1, nz, 'comp_i', [0, 1], 'Radius', 0.5, ...
                 'name', 'Injector', 'type', 'rate', 'Val', rate, 'sign', 1);
W(1).components = [0.0, 0.95, 0.05, 0.0, 0.0];  % H2-rich injection

%% Model Setup: Compositional Model with Bacterial Growth
arg = {G, rock, fluid, compFluid, 'water', false, 'oil', true, 'gas', true, ...
       'bacteriamodel', true, 'diffusioneffect', false, 'liquidPhase', 'O', 'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
model.outputFluxes = false;

%% Initial Conditions
% Temperature and initial saturations
T0 = 317.5;                % Initial temperature (K)
s0 = [0.8, 0.2];           % Initial saturations (Sw = 1)
z0 = [0.8, 0.0, 0.006, 0.018, 0.176];  % Initial composition: H2O, H2, CO2, N2, CH4
Phydro0=rhow*norm(gravity).*G.cells.centroids(:,3);
% Initialize state with bacterial concentration
nbact0 = 10^6;
 state0 = initCompositionalStateBacteria(model, Phydro0, T0, s0, z0, nbact0);
% state0.s =state0.s.*0+s0;

%% Time Stepping and Schedule
% Define schedule and solver
deltaT = T / niter;
schedule = simpleSchedule(repmat(deltaT, 1, niter), 'bc', [], 'src', [], 'W', W);
nls = NonLinearSolver('useRelaxation', true);

% Run simulation
[~, states, report] = simulateScheduleAD(state0, model, schedule);

%% Plotting Results
time = 0;
figure;
for i = 1:niter
    % Reshape spatial data for plotting
    x = G.cells.centroids(:, 1);
    z = G.cells.centroids(:, 3);
    X = reshape(x, [nx, nz]);
    Z = reshape(z, [nx, nz]);
    zH2 = reshape(states{i}.components(:, 2), [nx, nz]);
    zH2O = reshape(states{i}.components(:, 1), [nx, nz]);
    zCO2 = reshape(states{i}.components(:, 3), [nx, nz]);
    Sw = reshape(states{i}.s(:, 1), [nx, nz]);
    Sg = reshape(states{i}.s(:, 2), [nx, nz]);
    Pres = reshape(states{i}.pressure, [nx, nz]);
    nbacteria = reshape(states{i}.nbact, [nx, nz]);

    % Plot saturation, gas fraction, bacterial concentration, and H2 mole fraction
    subplot(2, 2, 1);
    contourf(X, Z, Sw, 60, 'EdgeColor', 'auto');
    colorbar; colormap('jet');
    clim([0 1]);
    title('Water Saturation');
    ylabel('Depth (m)');
    axis equal; axis([0 Lx depth_res depth_res + H]);
    set(gca, 'Ydir', 'reverse');

    subplot(2, 2, 2);
    contourf(X, Z, Sg, 60, 'EdgeColor', 'auto');
    colorbar; colormap('jet');
    clim([0 1]);
    title('Gas Saturation');
    ylabel('Depth (m)');
    axis equal; axis([0 Lx depth_res depth_res + H]);
    set(gca, 'Ydir', 'reverse');

    subplot(2, 2, 3);
    contourf(X, Z, nbacteria, 60, 'EdgeColor', 'auto');
    colorbar; colormap('jet');
    title('Bacterial Concentration');
    ylabel('Depth (m)');
    axis equal; axis([0 Lx depth_res depth_res + H]);
    set(gca, 'Ydir', 'reverse');

    subplot(2, 2, 4);
    contourf(X, Z, zH2, 60, 'EdgeColor', 'auto');
    colorbar; colormap('jet');
    clim([0 0.8]);
    title('H2 Mole Fraction');
    ylabel('Depth (m)');
    axis equal; axis([0 Lx depth_res depth_res + H]);
    set(gca, 'Ydir', 'reverse');

    time = time + deltaT;
%     suptitle(sprintf('Injection Duration = %.2f days', convertTo(time, day)));
    pause(0.001);
end
