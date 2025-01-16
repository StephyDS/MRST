classdef DecayBactRateSRC <  StateFunction
    % The bacterial decay rate, given per cell
    properties           
    end
    
    methods
        function gp = DecayBactRateSRC(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'Psi_decay';
        end

        function Psidecay = evaluateOnDomain(prop, model, state)
            Psidecay = 0;
            bbact = model.ReservoirModel.b_bact;
            namecp = model.ReservoirModel.getComponentNames();
            pv = model.ReservoirModel.operators.pv;  % Reservoir porosity
            %pv = model.ReservoirModel.rock.poro;
            % Find indices for 'H2'
            idx_H2 = find(strcmp(namecp, 'H2'), 1);     % Locate 'H2'
            if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase && (~isempty(idx_H2))                                    
                s = model.ReservoirModel.getProps(state, 's');                        
                nbact = model.ReservoirModel.getProps(state, 'nbact');
                L_ix = model.ReservoirModel.getLiquidIndex();
                sL = s{L_ix};
                Psidecay = pv.*bbact.*nbact.*(nbact.*sL);           
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
