classdef DispersMolecDiffusFlowDiscretization < FlowDiscretization
    % DispersMolecDiffusFlowDiscretization
    % Discretization and state function grouping for bio-chemical flow
    % within a compositional model with molecular diffusion.

    properties
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function props = DispersMolecDiffusFlowDiscretization(model)
            % Constructor: inherit base FlowDiscretization properties
            props = props@FlowDiscretization(model);
            useMolDiff = isprop(model, 'molecularDiffusion') && model.molecularDiffusion;
            useMolDisp = isprop(model, 'molecularDispersion') && model.molecularDispersion;
            if (useMolDiff && ~useMolDisp)
                props = props.setStateFunction('ComponentTotalFlux', ...
                    ComponentTotalFluxMolecularDiffusion(model));
                props = props.setStateFunction('ComponentPhaseMolecularDiffFlux', ...
                    ComponentPhaseMolecularDiffFlux(model));
                props = props.setStateFunction('ComponentTotalMolecularDiffFlux', ...
                    ComponentTotalMolecularDiffFlux(model));
            elseif (~useMolDiff && useMolDisp)
                props = props.setStateFunction('ComponentTotalFlux', ...
                    ComponentTotalFluxMolecularDispersion(model));
                props = props.setStateFunction('ComponentPhaseMolecularDispFlux', ...
                    ComponentPhaseMolecularDispFlux(model));
                props = props.setStateFunction('ComponentTotalMolecularDispFlux', ...
                    ComponentTotalMolecularDispFlux(model));
            elseif (useMolDiff && useMolDisp)
                props = props.setStateFunction('ComponentTotalFlux', ...
                    ComponentTotalFluxMolecularDispersionDiffusion(model));
                props = props.setStateFunction('ComponentPhaseMolecularDispFlux', ...
                    ComponentPhaseMolecularDispFlux(model));
                props = props.setStateFunction('ComponentTotalMolecularDispFlux', ...
                    ComponentTotalMolecularDispFlux(model));
                props = props.setStateFunction('ComponentPhaseMolecularDiffFlux', ...
                    ComponentPhaseMolecularDiffFlux(model));
                props = props.setStateFunction('ComponentTotalMolecularDiffFlux', ...
                    ComponentTotalMolecularDiffFlux(model));
            else
                props = props.setStateFunction('ComponentTotalFlux', ComponentTotalFlux(model));
            end


        end

        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Call parent method for standard component conservation
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
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