classdef GrowthBactRateSRC < StateFunction
    % Bacterial growth rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   gr = GrowthBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the bacterial growth rate in each grid cell based on:
    %   - H2 and CO2 concentrations (Monod kinetics)
    %   - Bacterial population density
    %   - Liquid phase saturation and density
    %   - Pore volume
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

    properties
        % No additional properties needed
    end

    methods
        function gp = GrowthBactRateSRC(model, varargin)
            % Constructor for bacterial growth rate calculator

            % Initialize base class
            gp@StateFunction(model, varargin{:});

            % Declare dependencies - kinetic term only (no mass, volume, saturation)
            gp = gp.dependsOn({'x'}, 'state');        % Mole fractions for Monod kinetics

            gp.label = '\Psi_{growth}'; % LaTeX-style label
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            % Compute specific growth rate coefficient (kinetic term only)
            %
            % Returns: Psigrowthmax * axH2 * axsub [1/s]
            % Used with BacterialMass: source = Psigrowth * BacterialMass
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct containing fields
            %
            % RETURNS:
            %   Psigrowth - Specific growth rate coefficient [1/s]

            % Initialize with zero growth rate
            Psigrowth = 0;

            % Get component names and indices
            rm = model.ReservoirModel;
            bcrm=rm.biochemFluid;
            namecp = rm.getComponentNames();
            idx_H2 = find(strcmpi(namecp, bcrm.rH2), 1);
            idx_sub = find(strcmpi(namecp, bcrm.rsub), 1);

            % Check if bacterial modeling is active and components exist
            if strcmp(bcrm.metabolicReaction, 'MethanogenicArchae')
                if ~(rm.bacteriamodel && rm.liquidPhase && ~isempty(idx_H2) && ~isempty(idx_sub))
                    return;
                end
            end

            % Get required properties
            x = rm.getProp(state, 'x');

            % Extract liquid phase component mole fractions
            if iscell(x)
                xH2 = x{idx_H2};
                xsub = x{idx_sub};
            else
                xH2 = x(:, idx_H2);
                xsub = x(:, idx_sub);
            end

            % Get growth parameters
            alphaH2 = bcrm.alphaH2;
            alphasub = bcrm.alphasub;
            Psigrowthmax = bcrm.Psigrowthmax;

            % Calculate Monod kinetics: growth rate coefficient (1/s)
            axH2 = xH2 ./ (alphaH2 + xH2);
            axsub = xsub ./ (alphasub + xsub);

            % Specific growth rate (kinetic term only; scaling by BacterialMass in flow)
            Psigrowth = Psigrowthmax .* axH2 .* axsub;
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