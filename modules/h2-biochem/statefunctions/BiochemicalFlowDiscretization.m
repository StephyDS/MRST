classdef BiochemicalFlowDiscretization < FlowDiscretization
    % BiochemicalFlowDiscretization
    % Discretization and state function grouping for bio-chemical flow
    % within a compositional model with microbial growth.

    properties
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function props = BiochemicalFlowDiscretization(model)
            % Constructor: inherit base FlowDiscretization properties
            props = props@FlowDiscretization(model);
            props = props.setStateFunction('ComponentTotalFlux', ...
                ComponentTotalFluxForBio(model));

            % Add dispersion state functions
            if (isprop(model, 'molecularDispersion') && model.molecularDispersion) || ...
                    (isprop(model, 'molecularDiffusion') && model.molecularDiffusion)
                props = props.setStateFunction('DispersiveDiffusivity', DispersiveDiffusivity(model));
                props = props.setStateFunction('DispersiveTransmissibility', ...
                    DispersiveTransmissibility(model));
                props = props.setStateFunction('ComponentPhaseDispFlux', ...
                    ComponentPhaseDispFlux(model));
                props = props.setStateFunction('ComponentTotalDispFlux', ...
                    ComponentTotalDispFlux(model));
            end

            if model.bacteriamodel
                % Ensure Porosity exists in FlowDiscretization for diffusion (also when no clogging)
                props = props.setStateFunction('Porosity', PorosityFromRock(model));

                % Set up transmissibility and porosity functions
                if model.dynamicFlowTrans
                    % Dynamic transmissibility computed each nonlinear iteration
                    props = props.setStateFunction('Permeability', BactPermeability(model));
                    props = props.setStateFunction('Transmissibility', ...
                        DynamicFlowTransmissibility(model, 'Permeability'));
                end

                if model.dynamicFlowPv
                    % Dynamic porosity / pore volume
                    props = props.setStateFunction('Porosity', BactPorosity(model));
                    props = props.setStateFunction('PoreVolume', ...
                        DynamicFlowPoreVolume(model, 'Porosity'));
                end

                if model.bactDiffusion
                    props = props.setStateFunction('MicrobialDiffusivity', MicrobialDiffusivity(model));
                    props = props.setStateFunction('MicrobialTransmissibility', ...
                        DynamicFlowTransmissibility(model, 'MicrobialDiffusivity'));
                    props = props.setStateFunction('BactFlux', DiffusiveBactFlux(model));
                end
            end
        end

        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Call parent method for standard component conservation
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
        end

        %-----------------------------------------------------------------%
        function [acc, bflux, name, type] = bacteriaConservationEquation(fd, model, state, state0, dt)
            % Computes bacterial mass conservation equation
            bactmass  = model.getProp(state, 'BacterialMass');
            bactmass0 = model.getProp(state0, 'BacterialMass');

            % Accumulation term: d(bactmass)/dt
            acc = (bactmass - bactmass0) ./ dt;

            % Add microbial diffusion contributions if present
            bflux=[];
            if (model.bactDiffusion)
                flowState = fd.buildFlowState(model, state, state0, dt);
                bflux = model.getProp(flowState, 'BactFlux');
            end

            % Convert to cell arrays for AD framework
            acc   = {acc};
            bflux = {bflux};

            % Output variable names and types
            name = 'bacteria';
            type = 'cell';
        end

        %-----------------------------------------------------------------%
        function c = getPopulationMass(model, state, extra)
            % Compute microbial mass in each phase
            c = component.getPhaseComposition(model, state);

            if nargin < 4
                rho = model.getProp(state, 'Density');
            else
                rho = extra.rho;
            end

            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = rho{ph} .* c{ph};
                end
            end
        end
    end
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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