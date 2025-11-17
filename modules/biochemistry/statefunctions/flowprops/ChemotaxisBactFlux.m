classdef ChemotaxisBactFlux < StateFunction   
    % bacterial Chemotaxis flux
    properties 

    end
    
    methods
        function gp = ChemotaxisBactFlux(model, varargin)                                  
            gp@StateFunction(model);
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            
            gp = gp.dependsOn({'s'}, 'state');
            gp = gp.dependsOn({'x'}, 'state');
            gp.label = 'Flux_{bioChemotaxis}';
        end
        function fluxchemobact = evaluateOnDomain(prop, model, state)
            % Get dependencies
            avg = model.operators.faceAvg;
            nbact = model.getProps(state, 'nbact');
            s = model.getProps(state, 's');                
            x = model.getProps(state, 'x');                
            L_ix = model.getLiquidIndex();
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');

            namecp =model.compFluid.names();
            indH2=find(strcmp(namecp, 'H2'));

            if model.liquidPhase && (~isempty(indH2))
                if iscell(s)     
                    sL = s{L_ix};  
                    rhoL = rho{L_ix};  
                else
                    sL = s(:, L_ix);
                    rhoL = rho(:, L_ix); 
               end
          
                if iscell(x)                  
                    xH2=x{indH2};
                else
                    xH2=x(:,indH2);
                end

                numChemBact=max(0.0,(xH2-model.xch_seuil));
                ChemBact=model.Xch_max.*numChemBact./(model.Kch+numChemBact);
                Diffch = avg(ChemBact.*sL.*pv.*rhoL.*nbact);
         
                fluxchemobact = -Diffch.*model.operators.Grad(xH2);

            else
                fluxchemobact =0;
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