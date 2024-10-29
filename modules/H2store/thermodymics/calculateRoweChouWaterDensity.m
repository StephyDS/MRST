function rho = calculateRoweChouWaterDensity(T, S)
% This function computes the water density using the Rowe correlation.
%
% SYNOPSIS:
%   rho = calculateRoweChouWaterDensity(T, S)
%
% DESCRIPTION:
%   This function calculates the water density based on the Rowe correlation,
%   incorporating temperature (T) and salinity (S). The correlation is derived
%   from empirical coefficients obtained for water properties.
%
% INPUTS:
%   T - Temperature in Kelvin (K)
%   S - Salinity in grams per liter (g/L)
%
% OUTPUTS:
%   rho - Water density in kilograms per cubic meter (kg/m^3)
%
% REFERENCE:
%   Safari Raad, Seyed Mostafa, Ranjbar, Ehsan, Hassanzadeh, Hassan,
%   Leonenko, Yuri, 2023. Hydrogen-brine mixture PVT data for reservoir
%   simulation of hydrogen storage in deep saline aquifers.
%   Int. J. Hydrog. Energy (ISSN: 0360-3199) 48 (2), 696–708.

% Coefficients
a1 = 5.916365 - 0.01035794 .* T + 0.9270048e-5 .* T.^2 - 1127.522 ./ T + 100674.1 ./ T.^2;
a2 = 0.520491e-2 - 0.10482101e-4 .* T + 0.8328532e-8 .* T.^2 - 1.1702939 ./T + 102.2783 ./T.^2;
a3 = 0.118547e-7 - 0.6599143e-10 .* T;
a4 = -2.5166 + 0.0111766 .* T - 0.170522e-4 .* T.^2;
a5 = 2.84851 - 0.0154305 .* T + 0.223982e-4 .* T.^2;
a6 = -0.0014814 + 0.829639e-5 .* T - 0.12469e-7 .* T.^2;
a7 = 0.0027141 - 0.15391e-4 .* T + 0.22655e-7 * T.^2;
a8 = 0.62158e-6 - 0.40075e-8 .* T + 0.65972e-11 .* T.^2;

% Calculate density using Rowe correlation
rho_inverse = (a1 - a2 * pi - a3 .* pi^2 + a4 * S + a5 .* S.^2 - a6 *pi.* S ...
    - a7 * pi * S.^2 - 0.5 .* a8 * pi^2 * S).*10^-3;

% Density is the inverse of rho_inverse
rho = 1 ./ rho_inverse; % Convert to kg/m^3
end
