classdef BiochemistryGenericFacilityModel < GenericFacilityModel
    % BiochemistryGenericFacilityModel
    % Generic facility model for biochemistry simulations, including
    % bacterial growth and decay.

    properties
        bacterialFormulation = 'bacterialmodel';  % Formulation for bacterial transport
    end

    methods
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, useDefaults)
            % Set up state function groupings using parent, and add
            % biochemistry-specific functions

            if nargin < 2
                useDefaults = isempty(model.FacilityFlowDiscretization);
            end

            % Base facility groupings
            model = setupStateFunctionGroupings@GenericFacilityModel(model, useDefaults);

            % Add biochemistry-specific state functions only when bacteriamodel is active
            rm = model.ReservoirModel;
            if ~isempty(rm) && isprop(rm, 'bacteriamodel') && rm.bacteriamodel
                ffd = model.FacilityFlowDiscretization;
                % Note: BacterialMass is already registered in ReservoirModel's PVTPropertyFunctions
                ffd = ffd.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                ffd = ffd.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                ffd = ffd.setStateFunction('BactConvRate', BactConvertionRate(model));
                model.FacilityFlowDiscretization = ffd;
            end
        end

        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            % Initialize AD state from double state
            state = initStateAD@GenericFacilityModel(model, state, vars, names, origin);
        end

        %-----------------------------------------------------------------%
        function names = getBasicPrimaryVariableNames(model)
            % Get names of primary variables
            names = getBasicPrimaryVariableNames@GenericFacilityModel(model);

            % If bacterial variables are disabled, return immediately
            if strcmpi(model.bacterialFormulation, 'none') || ...
                    strcmpi(model.primaryVariableSet, 'none')
                return
            end
        end

        %-----------------------------------------------------------------%
        function [variables, names, map] = getBasicPrimaryVariables(model, wellSol)
            % Return facility primary variables
            [variables, names, map] = getBasicPrimaryVariables@GenericFacilityModel(model, wellSol);
        end

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get field name and index for a given variable
            [fn, index] = getVariableField@GenericFacilityModel(model, name, varargin{:});
        end

        %-----------------------------------------------------------------%
        function src = getBacteriaSources(model, fd, state, state0, dt)
            % Growth/decay source: src = (g - d) * BacterialMass
            % where g, d are kinetic rates (1/s) applied to mass directly
            rm = model.ReservoirModel;
            if isempty(rm) || ~isprop(rm, 'bacteriamodel') || ~rm.bacteriamodel
                src = 0;
                return;
            end

            reg = 1.0e-10;
            flowState = fd.buildFlowState(model, state, state0, dt);
            psigrowth = model.getProps(flowState, 'PsiGrowthRate');  % Psigrowthmax * axH2 * axsub [1/s]
            psidecay  = model.getProps(flowState, 'PsiDecayRate');   % bbact * nbact [1/s]
            bmass     = rm.PVTPropertyFunctions.get(rm, state, 'BacterialMass');  % pv * S_l * rho_l * nbact [kg]

            % Direct (g-d)*mass formulation using BacterialMass
            src_growthdecay = (psigrowth - psidecay) .* bmass - reg .* bmass;

            % ===== NEW: Well bacteria source (advective transport) =====
            map   = model.getProp(state, 'FacilityWellMapping');
            rm = model.ReservoirModel;
            if ~isempty(map.cells)
                q_ph  = model.getProp(state, 'PhaseFlux');
                rho = rm.PVTPropertyFunctions.get(rm, state, 'Density');
                nbact = rm.getProp(state, 'nbact');
                L_ix  = rm.getLiquidIndex();

                % Liquid phase flux per perforation (positive = injection)
                q_l   = q_ph{L_ix};
                % Liquid density in perforated cells
                rho_l = rho{L_ix};
                % Bacteria mass flux: ρ_l * q_l * ω
                rho_perf = rho_l(map.cells);
                % Bacteria mass flux at each perforation
                q_bact = rho_perf .* q_l .* nbact(map.cells);
                % Injectors: no bacteria injected → set to 0
                q_bact(q_l > 0) = 0;

                % Sum perforation contributions to cells (producers give negative)
                % Sparse summation to cells (AD‑compatible)
                nc = rm.G.cells.num;
                S = sparse(map.cells, (1:numel(map.cells))', 1, nc, numel(map.cells));
                src_well = S * q_bact;
            else
                src_well = 0;
            end
            % ==============================================================

            src = src_growthdecay + src_well;
        end
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Return facility equations including parent contributions
            [eqs, names, types, state] = ...
                getModelEquations@GenericFacilityModel(model, state0, state, dt, drivingForces);
        end

        %-----------------------------------------------------------------%
        function [values, tolerances, names, evaluated] = getFacilityConvergenceValues(model, problem, varargin)
            % Get convergence values for facility
            [values, tolerances, names, evaluated] = ...
                getFacilityConvergenceValues@GenericFacilityModel(model, problem, varargin{:});
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