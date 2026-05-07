classdef ChemotaxisBactFlux < StateFunction
    % Bacterial chemotaxis flux computation for compositional simulations
    %
    % SYNOPSIS:
    %   Jchemotaxis = ChemotaxisBactFlux(model, 'property1', value1, ...)
    %
    % Jchemotaxis = Kchemo * nbact *(1 - nbact) * Grad(xH2)
    % (volume-filling effect, generalized Keller-Segel model)
    %
    % DESCRIPTION:
    %   Computes the bacterial Chemotaxis flux in each grid cell, accounting for:
    %   - Bacterial concentration (nbact)
    %   - Liquid phase saturation
    %   - Pore volume
    %   - Phase densities
    %   - Presence of required components H2
    %   - displacement along with the H2 concentration gradient
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with bacterial modeling enabled
    %
    % OPTIONAL PARAMETERS:
    %   None
    %
    % RETURNS:
    %   Class instance ready for use in simulation
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional
    %
    %Sources: https://www.cirm-math.fr/RepOrga/2260/Slides/cours_Ribot.pdf

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
            op=model.operators;
            bcrm=model.biochemFluid;
            nbact = model.getProps(state, 'nbact');
            s = model.getProps(state, 's');
            x = model.getProps(state, 'x');
            L_ix = model.getLiquidIndex();
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');

            namecp = model.getComponentNames();
            nbioreact=bcrm.nbioreact;
            % Initialize with zero growth rate
            fluxchemobact=cell(1,nbioreact);
            [fluxchemobact{:}] = deal(0);

            % Check if bacterial modeling is active and components exist
            if ~(model.bacteriamodel && model.liquidPhase)
                return;
            end


            for i=1:nbioreact
                indH2= find(strcmpi(namecp, bcrm.rH2(i)), 1);

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
                    if iscell(nbact)
                        nbacti=nbact{i};
                    else
                        nbacti=nbact(:,i);
                    end

                    % Calculate effective volume with safeguards
                    if iscell(sL)
                        Voln = max(sL{1}, 1.0e-8) .* rhoL{1};
                    else
                        Voln = max(sL, 1.0e-8) .* rhoL;
                    end
                    Voln = max(Voln, 1.0e-8);


                    %deltaChiBact=max(0.0, xH2-bcrm.xch_seuil(i));
                    %ChiBact=bcrm.Xch_max(i).*deltaChiBact ./ (bcrm.Kch(i)+deltaChiBact);
                    ChiBact=bcrm.xch_seuil(i);
                    Diffch = op.faceAvg(ChiBact.*pv.*Voln.*nbacti.*(1.0-nbacti));
                    fluxchemobact{i} = Diffch.*model.operators.Grad(xH2);            
                end
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
