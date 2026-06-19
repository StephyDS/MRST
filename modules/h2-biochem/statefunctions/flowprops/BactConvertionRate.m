classdef BactConvertionRate < StateFunction
    % Bacterial conversion rate model for compositional simulations
    %
    % SYNOPSIS:
    %   gf = BactConvertionRate(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes component conversion rates due to bacterial growth.
    %   Uses kinetic growth rate (without mass) multiplied by current bacterial mass:
    %   qbiot = (growth_rate * nbact) * stoichiometric_coefficient
    %   This ensures consistent treatment of (g-d)*mass formulation in the flow.
    %
    % SEE ALSO:
    %   CompositionalModel, EquationsCompositional

    properties
        % No additional properties needed - all parameters come from model
    end

    methods
        function gp = BactConvertionRate(model, varargin)
            % Constructor for bacterial conversion rate calculator
            gp@StateFunction(model, varargin{:});

            % Define dependencies - computes kinetics directly, no pre-computed rates
            gp = gp.dependsOn({'Density'}, 'PVTPropertyFunctions');

            % Set label for output
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            % Compute bacterial conversion rate for each component
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Compositional model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   qbiot - Cell array of conversion rates per component [kg/s]

            % Initialize with zeros
            ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents();
            qbiot = cell(ncomp, 1);
            [qbiot{:}] = deal(0);

            % Check if bacterial modeling is active
            rm = model.ReservoirModel;
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end
            bcrm=rm.biochemFluid;

            % Get component indices
            compNames = rm.EOSModel.getComponentNames();
            if strcmp(bcrm.metabolicReaction, 'MethanogenicArchae')
                idxH2  = find(strcmpi(compNames, bcrm.rH2),1);
                idxsub = find(strcmpi(compNames, bcrm.rsub),1);

                % Validate required components
                if isempty(idxH2) || isempty(idxsub)
                    warning('Bacterial model requires H2 and CO2 components');
                    return;
                end
            end

            try
                % Get required state variables
                rho = rm.PVTPropertyFunctions.get(model.ReservoirModel, state, 'Density');
                L_ix = rm.getLiquidIndex();
                x = rm.getProp(state, 'x');

                % Extract liquid phase properties
                if iscell(x)
                    xH2 = x{idxH2};
                    xsub = x{idxsub};
                    rhoL = rho{L_ix};
                else
                    xH2 = x(:, idxH2);
                    xsub = x(:, idxsub);
                    rhoL = rho(:, L_ix);
                end

                % Ensure non-zero liquid saturation
                rhoL = max(value(rhoL), 1.0e-8);

                % Get model parameters
                alphaH2 = bcrm.alphaH2;
                alphasub = bcrm.alphasub;
                Psigrowthmax = bcrm.Psigrowthmax;
                Y_H2 = bcrm.Y_H2;
                gammak = rm.gammak;
                nbactMax = bcrm.nbactMax;
                mc = rm.EOSModel.CompositionalMixture.molarMass;

                % Calculate Monod kinetics (specific growth rate)
                axH2 = xH2 ./ (alphaH2 + xH2);
                axsub = xsub ./ (alphasub + xsub);
                growth_kinetic = Psigrowthmax .* axH2 .* axsub;

                % Theory: q_i = Φ * γ_i^H2 * ψ_growth * S_l / Y_H2
                % where γ_i^H2 = n0 * γ_i * M_i / γ_H2
                % and ψ_growth = growth_kinetic * nbact (biomass growth rate)

                % Compute stoichiometric factor for each component
                % γ_i^H2 = nbactMax * γ_i * M_i / γ_H2
                stoich_factor = nbactMax .* mc ./ abs(gammak(idxH2));

                % Component source: q_i = pv * γ_i^H2 * (growth_kinetic * nbact) * S_l / Y_H2
                % Rewrite: q_i = pv * S_l * nbact / Y_H2 * (γ_i^H2 * growth_kinetic)
                qbiot_temp =  growth_kinetic ./ (Y_H2 .*rhoL);

                for c = 1:ncomp
                    qbiot{c} = gammak(c) .* stoich_factor(c) .* qbiot_temp;
                end

            catch ME
                warning('Bacterial conversion rate calculation failed: %s', ME.message);
                [qbiot{:}] = deal(0);
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