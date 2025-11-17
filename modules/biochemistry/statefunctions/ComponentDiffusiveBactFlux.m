classdef ComponentDiffusiveBactFlux < StateFunction  
    % bacterial diffusive flux
    properties
    end

    methods
        function gp = ComponentDiffusiveBactFlux(model, varargin)
            gp@StateFunction(model);
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            
            gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'Flux_{bio}';
        end

        function fluxbact = evaluateOnDomain(prop, model, state)
            % Get dependencies
            avg = model.operators.faceAvg;
            nbact = model.getProps(state, 'nbact');
            s = model.getProps(state, 's');
            L_ix = model.getLiquidIndex();
            rho = prop.getEvaluatedExternals(model, state, 'Density');
            poro= model.rock.poro;
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
           
           fluxbact = model.AutoDiffBackend.initVariablesAD(zeros(model.G.cells.num, 1)); % AD-safe init
           

            if model.bacteriamodel&& model.liquidPhase
                if iscell(s)
                    sL = s{L_ix};
                    rhoL = rho{L_ix};
                else
                    sL = s(:, L_ix);
                    rhoL = rho(:, L_ix);
                end
                tau_mq=(sL.*poro).^(7/3).*poro.^(-2);
                Diffb = avg(model.bact_diff.*tau_mq.*sL.*pv.*rhoL);

                fluxbact = -Diffb.*model.operators.Grad(nbact);
            end
        end
    end
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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