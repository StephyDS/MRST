classdef BactConvertionRate < StateFunction
    % Bacterial conversion rate model for compositional simulations.
    %
    % SYNOPSIS:
    %   bcr = BactConvertionRate(model)
    %   bcr = BactConvertionRate(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Computes component conversion rates due to bacterial growth for
    %   methanogenic archaea. The conversion rate follows the stoichiometry
    %   of the metabolic reaction:
    %
    %       4 H2 + CO2 -> CH4 + 2 H2O
    %
    %   The conversion rate for each component c is computed as:
    %
    %       qbiot(c) = γ_c * (bmass / (ρL * Y_H2))
    %
    %   where:
    %       γ_c   = Stoichiometric coefficient for component c
    %       bmass = Bacterial mass concentration [kg/m^3]
    %       ρL    = Liquid phase density [kg/m^3]
    %       Y_H2  = Yield coefficient for H2 [kg_biomass/kg_H2]
    %
    %   The stoichiometric coefficients are normalized by the H2 coefficient
    %   and scaled by the bacterial carrying capacity (nbactMax).
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance with:
    %           .bacteriamodel (true)
    %           .biochemFluid.metabolicReaction ('MethanogenicArchae')
    %           .biochemFluid.rH2 (reactant component name, e.g., 'H2')
    %           .biochemFluid.rsub (substrate component name, e.g., 'CO2')
    %           .biochemFluid.Y_H2 (yield coefficient)
    %           .biochemFluid.nbactMax (carrying capacity)
    %           .gammak (stoichiometric coefficients)
    %
    % RETURNS:
    %   qbiot - Cell array of conversion rates for each component [kg/s]
    %
    % SEE ALSO:
    %   StateFunction, BiochemistryModel, ComponentSource

    properties
        % No additional properties needed - all parameters come from model
    end

    methods
        %-----------------------------------------------------------------%
        function bcr = BactConvertionRate(model, varargin)
            % Constructor for bacterial conversion rate calculator.
            %
            % PARAMETERS:
            %   model - Compositional model with bacterial physics
            %   varargin - Optional property/value pairs
            %
            % RETURNS:
            %   bcr - BactConvertionRate instance

            bcr@StateFunction(model, varargin{:});

            % Dependencies: pore volume, saturation, and bacterial concentration
            % Note: Direct calculation avoids redundant Density computation through BacterialMass
            bcr = bcr.dependsOn('BacterialMass', 'PVTPropertyFunctions');
            bcr = bcr.dependsOn('PsiGrowthRate', 'state');

            bcr.label = 'Q_biot';
        end

        %-----------------------------------------------------------------%
        function qbiot = evaluateOnDomain(bcr, model, state)
            % Compute bacterial conversion rate for each component.
            %
            % PARAMETERS:
            %   model - Compositional model instance
            %   state - Reservoir state structure
            %
            % RETURNS:
            %   qbiot - Cell array qbiot{c} containing conversion rates
            %           for component c (ncells × 1), units [kg/s]

            % Get reservoir model and number of components
            rm = model.ReservoirModel;
            bcrm=rm.biochemFluid;
            ncomp = rm.EOSModel.getNumberOfComponents();
            nbioreact=bcrm.nbioreact;
            qbiot = cell(ncomp, 1);
            [qbiot{:}] = deal(0);

            % Early return if bacterial model is not active
            if ~(rm.bacteriamodel && rm.liquidPhase)
                return;
            end

            % Validate metabolic reaction type
            if ~strcmp(bcrm.metabolicReaction, 'MethanogenicArchae')
                error('BactConvertionRate:UnsupportedReaction', ...
                    'Unsupported metabolic reaction: %s', bcrm.metabolicReaction);
            end

            % Get component names and find indices of H2 and CO2
            compNames = rm.EOSModel.getComponentNames();
            for i=1:nbioreact
                if strcmp(bcrm.metabolicReaction(i), 'MethanogenicArchae') || ...
                        strcmp(bcrm.metabolicReaction(i), 'AcetogenicBacteria')
                    idxH2  = find(strcmpi(compNames, bcrm.rH2(i)),1);
                    idxsub = find(strcmpi(compNames, bcrm.rsub(i)),1);

                    % Validate required components exist
                    if isempty(idxH2)
                        error('BactConvertionRate:MissingH2', ...
                            'Required component H2 (%s) not found in the system', bcrm.rH2);
                    end
                    if isempty(idxsub)
                        error('BactConvertionRate:MissingCO2', ...
                            'Required component CO2 (%s) not found in the system', bcrm.rsub);
                    end
                end
            end

            % Get primary state variables directly (avoids redundant Density calculations)
            bmass = rm.PVTPropertyFunctions.get(rm, state, 'BacterialMass');
            psigrowth = model.getProps(state, 'PsiGrowthRate');

            % Get model parameters
            for i=1:nbioreact
                idxH2  = find(strcmpi(compNames, bcrm.rH2(i)),1);
                Y_H2 = bcrm.Y_H2(i);                      % Yield coefficient [kg_biomass/kg_H2]
                gamma = rm.gammak(i,:);                     % Stoichiometric coefficients
                nbactMax = bcrm.nbactMax(i);              % Carrying capacity [cells/m^3]
                mc = rm.EOSModel.CompositionalMixture.molarMass;  % Component molar masses [kg/mol]

                % Normalize stoichiometric coefficients by H2 coefficient
                % γ_c^H2 = (γ_c * m_c) / |γ_H2|
                % Scaled by nbactMax for consistency with bacterial mass units
                gamma_norm = nbactMax .* gamma .* mc ./ abs(gamma(idxH2));

                % Calculate base conversion rate
                % Direct calculation: qbiot = (pv * S_l * nbact) / Y_H2
                % Note: Density cancels in (pv * S_l * rho_l * nbact) / (rho_l * Y_H2)
                qbase = psigrowth{i}.*bmass{i} ./ Y_H2;

                % Apply normalized stoichiometric coefficients to each component
                for c = 1:ncomp
                    qbiot{c} =qbiot{c}+ gamma_norm(c) .* qbase;
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