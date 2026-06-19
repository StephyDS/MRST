classdef DecayBactRateSRC < StateFunction
    % Bacterial decay rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   decay = DecayBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the bacterial decay rate in each grid cell, accounting for:
    %   - Bacterial concentration (nbact)
    %   - Liquid phase saturation
    %   - Pore volume
    %   - Phase densities
    %   - Presence of required components (H2 and CO2)
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
        function gp = DecayBactRateSRC(model, varargin)
            % Constructor for bacterial decay rate calculator
            gp@StateFunction(model, varargin{:});

            % Define dependencies - kinetic coefficient based on current concentration
            gp = gp.dependsOn({'nbact'}, 'state');          % Bacterial concentration

            % Set label for output
            gp.label = 'Psi_{decay}';
        end

        function Psidecay = evaluateOnDomain(prop, model, state)
            % Compute density-dependent decay rate coefficient
            %
            % For density-dependent decay: b*nbact^2 = (b*nbact) * nbact
            % Returns: bbact * nbact [1/s]
            % Used with BacterialMass: source = Psidecay * BacterialMass / nbact^(scaling)
            % Or more directly: decay_contribution = Psidecay * nbact (when multiplied by mass)
            %
            % PARAMETERS:
            %   prop  - Property function instance
            %   model - Reservoir model instance
            %   state - State struct with fields
            %
            % RETURNS:
            %   Psidecay - Decay rate coefficient (depends on current nbact) [1/s]

            % Initialize with zeros
            Psidecay = 0;

            % Get model parameters
            rm = model.ReservoirModel;
            bcrm=rm.biochemFluid;
            bbact = bcrm.bbact;

            % Check if bacterial modeling is active
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end

            % Get component names and indices
            namecp = rm.getComponentNames();
            idx_H2 = find(strcmpi(namecp, bcrm.rH2), 1);
            idx_sub = find(strcmpi(namecp, bcrm.rsub), 1);

            % Validate required components
            if strcmp(bcrm.metabolicReaction, 'MethanogenicArchae')
                if isempty(idx_H2) || isempty(idx_sub)
                    return;
                end
            end

            % Get required state variables
            nbact = rm.getProp(state, 'nbact');

            % Compute density-dependent decay coefficient: b * nbact [1/s]
            % When multiplied by BacterialMass and divided by nbact, gives b*nbact^2
            nbact_pos = max(nbact, 0);
            Psidecay = bbact .* nbact_pos;
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