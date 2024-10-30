% MRST Example: Setting Up H₂-Brine Fluid Properties for Black-Oil Simulation
%
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
% Parameters for Temperature, Pressure, and Saturation Range
min_temp = 40;                % Minimum temperature in Celsius
max_temp = 60;                % Maximum temperature in Celsius
min_pressure = 1 * atm;       % Minimum pressure in Pa
max_pressure = 20 * atm;      % Maximum pressure in Pa
nbp = 10;                     % Number of pressure points
nbt = 10;                     % Number of temperature points
ms = 2;                       % Scaling factor for solubility
% Get the current script's directory
currentDir = fileparts(mfilename('fullpath'));
outputDisplay = false;
% Define the target output directory relative to the current directory
output_dir = fullfile(currentDir, 'data', 'UHS_BENCHMARK_RSRV');

% Check if the directory exists; if not, create it
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
% Generate H2O Component Table
comp_name = 'H2O';
disp(['Generating component table for: ', comp_name]);
[tab_H2O, status_H2O, file_path_H2O] = generateComponentTable(min_temp, max_temp, min_pressure, max_pressure, nbp, nbt, comp_name, output_dir, 'outputDisplay', false);

% Generate H2 Component Table
pause(0.5);  % Ensure smooth execution between commands
comp_name = 'H2';
disp(['Generating component table for: ', comp_name]);
[tab_H2, status_H2, file_path_H2] = generateComponentTable(min_temp, max_temp, min_pressure, max_pressure, nbp, nbt, comp_name, output_dir,'outputDisplay', false);

% Generate Solubility Table
disp('Generating solubility table...');
[tab_sol, status_sol, file_path_sol] = generateSolubilityTable(min_temp, max_temp, min_pressure, max_pressure, nbp, nbt, ms, output_dir);

% Configure and Write Fluid Properties (PVT) Tables
getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol, ...
    'rs', true, 'rv', true, ...
    'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE', 'PVTGFile', 'PVTGH2BRINE', ...
    'dir', 'UHS_BENCHMARK_RSRV');

disp('Writing fluid properties for miscible case with digas and evapoil...');
getFluidH2BrineProps(tab_H2O, tab_H2, tab_sol, 'rs', true, 'rv', true, 'PVTGFile', 'PVTGH2BRINE', 'PVTOFile', 'PVTOH2BRINE','PVDOFile', 'PVDOH2BRINE', 'dir', output_dir)
%% This 2ndpar generates gas-oil flow properties specific to a hydrogen-brine system, modeled over
%% three distinct regions. We exemplify the system with three rock types: caprock, bedrock, and storage rock (aquifer).
%% The output includes key properties: gas relative permeability (krG), oil relative permeability (krO),
%% and capillary pressure at the gas-oil contact (pcOG). 
%
%% Relative permeability is modeled quadratically to capture the flow dynamics between gas and oil phases
%% within each rock type, providing a more realistic representation of gas-brine interactions under 
%% varying reservoir conditions.
%% Generate SGOF (Gas-Oil Flow) Properties Table
disp('Generating gas-oil flow properties for hydrogen-brine system with three regions');
getFluidH2BrineSGOF('n', 100, 'plot', true, 'dir', output_dir, 'filename', 'SGOF_H2STORAGE.txt', ...
                      'units', 'metric', 'nreg', 3, 'sw_imm', [0.1, 0.1, 0.1], ...
                      'sg_imm', [0.1, 0.1, 0.1], 'c_a1', 2, 'c_a2', 2, 'c_a3', 1.5, ...
                      'Pe', [0.4, 0.4, 0.4], 'P_c_max', [9.5e4, 9.5e4, 9.5e4]);

% Display completion message
disp('H₂-Brine fluid properties setup for black-oil simulation is complete.');