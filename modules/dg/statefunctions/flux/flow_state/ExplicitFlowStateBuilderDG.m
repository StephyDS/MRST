classdef ExplicitFlowStateBuilderDG < ExplicitFlowStateBuilder & FlowStateBuilderDG
    methods
        function flowState = build(builder, fd, model, state, state0, dt, type)
            % Hybridize state. The base state is the implicit. Other
            % functions are then assigned.
            flowState = build@FlowStateBuilderDG(builder, fd, model, state, state0, dt, type);
            builder.explicitProps = {'s'};
            switch type
                case 'face'
                    flowState0 = state0.faceStateDG;
                case 'cell'
                    flowState0 = state0.cellStateDG;
                    builder.explicitFlux = {};
                    builder.explicitFlow = {'ComponentMobility', 'Mobility', 'GravityPermeabilityGradient'};
            end
            flowState = build@ExplicitFlowStateBuilder(builder, fd, model, flowState, flowState0, dt);
        end
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
