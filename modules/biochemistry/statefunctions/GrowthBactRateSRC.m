classdef GrowthBactRateSRC <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = GrowthBactRateSRC(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'x'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp = gp.dependsOn({'nbact'}, 'state');
            gp.label = 'Psi_growth';
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            Psigrowth = 0;
            namecp = model.ReservoirModel.getComponentNames();
            % Find indices for 'H2' and 'CO2'
            idx_H2 = find(strcmp(namecp, 'H2'));     % Locate 'H2'
            idx_CO2 = find(strcmp(namecp, 'CO2'));   % Locate 'CO2'

         if model.ReservoirModel.bacteriamodel&& model.ReservoirModel.liquidPhase && (~isempty(idx_H2)) && (~isempty(idx_CO2))
                          
             pv = model.ReservoirModel.rock.poro;
             x = model.ReservoirModel.getProps(state, 'x');           
             if iscell(x)                  
                 xH2 = x{idx_H2};     % Mole fraction of H2
                 xCO2 = x{idx_CO2};
             else
                 xH2 = x(:, idx_H2);  % Column corresponding to H2 in matrix
                 xCO2 = x(:, idx_CO2);  % Column corresponding to CO2 in matrix
             end

             s = model.ReservoirModel.getProps(state, 's');            
             nbact = model.ReservoirModel.getProps(state, 'nbact');

             L_ix = model.ReservoirModel.getLiquidIndex();
             sL = s{L_ix};
             alphaH2 = model.ReservoirModel.alphaH2;
             alphaCO2 = model.ReservoirModel.alphaCO2;
             Psigrowthmax = model.ReservoirModel.Psigrowthmax;
             % Calculate Psigrowth using H2 and CO2 mole fractions
             Psigrowth = pv.*Psigrowthmax .* (xH2 ./ (alphaH2 + xCO2)) ...
                        .* (xCO2 ./ (alphaCO2 + xH2)).*nbact.*sL;         
         end
        end
    end
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
