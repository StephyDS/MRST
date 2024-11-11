classdef BacterialMass < StateFunction & ComponentProperty
    % Mass of each component, in each phase.
    properties
    end
    
    methods
        function gp = BacterialMass(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'M_{bio}';
        end
        function mb = evaluateOnDomain(prop, model, state)
                 
            s = model.getProps(state, 's');
            nbact = model.getProps(state, 'nbact');

            pv = model.operators.pv;
            pv = model.rock.poro;
            L_ix = model.getLiquidIndex();
                
            if iscell(s)                  
                 sL = s{L_ix};    
             else
                 sL = s(:, L_ix);  
             end
            mb = pv.*nbact.*sL;
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
