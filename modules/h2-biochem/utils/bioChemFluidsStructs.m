function biochemfluids = bioChemFluidsStructs()
% Parameters for bioreactions
% 'metabolicReaction': name of bioreactions :'MethanogenicArchae','SulfateReducingBacteria'
%  Names of bioreactions reactants: 'rH2','rsub',
%  Names of bioreactions products: 'pH2O','p2',
%  Associated Stoichiometric coefficients: 'gamrH2','gamrsub','gampH2O','gamp2',
%  Yield coefficient:  'Y_H2',  1/mole(H2)
%  Half-saturation constant for H2: 'alphaH2', mol/mol
%  Half-saturation constant for substrate (CO2 or SO4); 'alphasub', mol/mol
%  Maximum growth rate: 'Psigrowthmax', 1/s
%  Constant decay rate:   'b_bact', 1/s
%  Maximum number of microorganisms per volume: 'nbactMax', 1/m3


biochemfluids = [...
    struct('metabolicReaction', 'MethanogenicArchae','rH2','H2','rsub','CO2',...
    'pH2O','H2O','p2','C1','gamrH2',-4,'gamrsub',-1,'gampH2O',2,'gamp2',1,...
    'Y_H2', 3.90875e11,'alphaH2',3.6e-7,'alphasub',1.98e-6,'Psigrowthmax',1.338e-4,...
    'bbact',2.35148e-6,'nbactMax',1e9,'bactdiff', 1.e-12, 'Xch_max', 1.e-9, ...
    'Kch', 5.e-6, 'xch_seuil', 1.e-7), ...
     ];
end
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
