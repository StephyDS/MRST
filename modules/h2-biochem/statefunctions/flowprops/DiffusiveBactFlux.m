classdef DiffusiveBactFlux < StateFunction
    % Bacterial diffusive flux computation for compositional simulations
    %
    % SYNOPSIS:
    %   jbactdiff = DiffusiveBactFlux(model, 'property1', value1, ...)
    %
    %Jbactdiff = - Kdiffbact * Grad(nbact)
    %
    % DESCRIPTION:
    %   Computes the bacterial Chemotaxis flux in each grid cell, accounting for:
    %   - Bacterial concentration (nbact)
    %   - Liquid phase saturation
    %   - Pore volume
    %   - Phase densities
    %   - Presence of required components H2
    %   - displacement along with the bacterial concentration gradient
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
    %Sources:
    % 1/ Phd thesis of Noura Eddaoui, Patterns in bioreactive transport 
    % in underground storage of hydrogen : impact of natural and induced
    % heterogeneity
    % 2/ Underground hydrogen storage: A comprehensive review, Davood
    % Zivar,https://doi.org/10.1016/j.ijhydene.2020.08.138
    %1/Impact of hydrogen on the hydrodynamic and bio-chemical
    % behavior, Hagemann & al., DOI 10.1007/s10596-015-9515-6



    properties
    end

    methods
        function gp = DiffusiveBactFlux(model, varargin)
            gp@StateFunction(model);
            gp = gp.dependsOn({'nbact'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');

            gp = gp.dependsOn({'s'}, 'state');
            gp.label = 'Flux_{bio}';
        end
        function fluxbact = evaluateOnDomain(prop, model, state)
            % Get dependencies
            op=model.operators;
            bcrm=model.biochemFluid;
            nbioreact=bcrm.nbioreact;
            nbact = model.getProps(state, 'nbact');
            s = model.getProps(state, 's');
            L_ix = model.getLiquidIndex();
            pv = model.PVTPropertyFunctions.get(model, state, 'PoreVolume');
            rho = model.PVTPropertyFunctions.get(model, state, 'Density');

           % Initialize with zero growth rate
            fluxbact=cell(1,nbioreact);
            [fluxbact{:}] = deal(0);
 
           % Check if bacterial modeling is active and components exist
            if ~(model.bacteriamodel && model.liquidPhase)
                return;
            end
            
            for i=1:nbioreact       
                if iscell(s)
                    sL = s{L_ix};
                    rhoL = rho{L_ix};
                else
                    sL = s(:, L_ix);
                    rhoL = rho(:, L_ix);
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
                
                Diffb = op.faceAvg(bcrm.bactdiff(i).*pv.*Voln);

                fluxbact{i} = -Diffb.*op.Grad(nbacti);
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
