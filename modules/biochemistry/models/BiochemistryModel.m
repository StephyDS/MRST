classdef BiochemistryModel <  GenericOverallCompositionModel 
%Bio-chemical model for compositional mixture with Hydrogen (H2)
%
% SYNOPSIS:
%   model = BioCompositionalModel(G, rock, fluid)
%   model = BioCompositionalModel(G, rock, fluid, compFluid)
%   model = BioCompositionalModel(..., 'pn1', vn1, ...)
%
% DESCRIPTION:
%   This model forms the basis for simulation of bio-chemical systems within compositional 
%   models. The model couple compositional modelto  bio-chemical reactions model. We use only 
%   microbial growth and decay of MOnod type.
%
% REQUIRED PARAMETERS:
%   G         - Simulation grid.
%
%   rock      - Valid rock used for the model.
%
%   fluid     - Fluid model used for the model.
%
%   compFluid - Compositional fluid mixture. Optional, defaults to H2O
%               mixture if left empty.
%
% OPTIONAL PARAMETERS:
%   'property' - Set property to the specified value.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   `ReservoirModel`, `ThreePhaseCompotitionalModel`
%

    properties

        % Boolean indicating if we are considering bio-chemical
        % effects (for debugging). Default = true.
        bactrial = true; 
        % Choice of bactrial primary variable. Possible values are 'dnbact'
        % for Plankton  and 'anbact' for biofilm. Default is 'nbact'.
        bacterialFormulation = 'bacterialmodel';

        % Compositional fluid mixutre. Default is a `Methanogenesis`
        % with H2-CO2-H2O-CH4
        compFluid 
        % Physical quantities and bounds
        Y_H2 = 3.9e11;  % Conversion factor for hydrogen consumption (moles/volume)
        
        %===========SDS MODIF===============================================
        %gammak = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0,0.0, -4];  % Stoichiometric coefficients: {'H2O'}    {'C1'}    {'N2'}    {'CO2'}    {'C2'}    {'C3'}    {'NC4'}    {'H2'}
        gammak = [];  %SDS MODIF 
        mol_diff = []; %zeros(8,2);
        
        % mol_diff = [ ...
        %               2.3e-9, 1.5e-5;  % Water (H2O)  
        %               2.6e-9, 1.6e-5;  % Methane (C1)
        %               2.1e-9, 1.8e-5;  % Nitrogen (N2)
        %               1.9e-9, 1.4e-5;  % Carbon Dioxide (CO2)                      
        %               3.2e-9, 2.5e-5;  % Ethane (C2)
        %               2.8e-9, 2.2e-5;  % Propane (C3)
        %               2.4e-9, 1.9e-5;  % Normal Butane (NC4)                      
        %               4.5e-9, 6.1e-5;  % Hydrogen (H2)
        %               ];
%===========SDS MODIF===============================================
        
        alphaH2 = 3.6e-7;
        alphaCO2 = 1.98e-6;
        Psigrowthmax = 1.338e-4;
        b_bact =2.3E-6;
        Db = 10^(-8)*meter/second
        bDiffusionEffect = false;
        moleculardiffusion = true;
        bacteriamodel = true;
        metabolicReaction='MethanogenicArchae' %SDS modif

    end
    
    methods
        %-----------------------------------------------------------------%
        function model = BiochemistryModel(G, rock, fluid, compFluid, varargin)
        % Class constructor. Required arguments are G, rock and fluid.
            model = model@GenericOverallCompositionModel(G, rock, fluid, compFluid);
            model = merge_options(model, varargin{:});
            % Check if compositional fluid is provided
            if nargin > 3 && isa(varargin{4}, 'CompositionalMixture')
                compFluid = varargin{4};
            end
            % Call parent constructor, but whithout rock because we will
            % set up operators later
            model = model.setupOperators();
            % Check if we have a water phase
            model.gas = true;
            if isfield(model.fluid, 'rhoO')
                model.oil = true;
            end

            %==================SDS MODIF=============================================
            namecp = compFluid.names();

            if model.moleculardiffusion
                indices = struct('H2', find(strcmp(namecp, 'H2')), ... 
                    'C1', find(strcmp(namecp, 'C1')), ... 
                    'CO2', find(strcmp(namecp, 'CO2')), ...
                    'H2O', find(strcmp(namecp, 'H2O')), ...
                    'N2', find(strcmp(namecp, 'N2')), ...
                    'C2', find(strcmp(namecp, 'C2')), ...
                    'C3', find(strcmp(namecp, 'C3')), ...
                    'NC4', find(strcmp(namecp, 'NC4')));
                
                fields = fieldnames(indices);
                nfields=numel(fields);
                model.mol_diff=zeros(nfields,2);

                coeffs = struct(...
                    'H2',  [4.5e-9, 6.1e-5], ...
                    'C1', [2.6e-9, 1.6e-5], ...
                    'H2O', [2.3e-9, 1.5e-5], ...
                    'CO2', [1.9e-9, 1.4e-5], ... 
                    'N2',  [2.1e-9, 1.8e-5], ... 
                    'C2',  [3.2e-9, 2.5e-5], ... 
                    'C3',  [2.8e-9, 2.2e-5], ... 
                    'NC4', [2.4e-9, 1.9e-5]);
                
                for i = 1:nfields
                    comp = fields{i};
                    indComp = indices.(comp);
                    
                    if ~isempty(indComp) && isfield(coeffs, comp)
                      model.mol_diff(indices.(comp),:)= coeffs.(comp);
                    end
                 end
            end

            % Set compositinal fluid
            if isempty(compFluid)
            
                % Default is Methanogenesis
                %compFluid = TableCompositionalMixture({'Hydrogen', 'Water','Nitrogen', 'CarbonDioxide', 'Methane'}, ...
                %   {'H2', 'Water', 'N2', 'CO2', 'CH4'});
                if strcmp(model.metabolicReaction,'MethanogenicArchae')
                    % Default is Methanogenesis
                    compFluid = TableCompositionalMixture({'Hydrogen', 'Water','Nitrogen', 'CarbonDioxide', 'Methane'}, ...
                    {'H2', 'Water', 'N2', 'CO2', 'C1'});
                    model.gammak = [-1.0, 0.5, 0.0, -0.25, 0.25];  %SDS MODIF 
               else
                    %send an error message: TO DO
               end
                %==================SDS MODIF=========================
                 %==================SDS MODIF=========================
            else
                %namecp = compFluid.names();
                ncomp=compFluid.getNumberOfComponents();
                model.gammak=zeros(1,ncomp);
                if strcmp(model.metabolicReaction,'MethanogenicArchae')
                    indH2=find(strcmp(namecp,'H2'));
                    indH2O= find(strcmp(namecp,'H2O'));
                    indCO2=find(strcmp(namecp,'CO2'));
                    indC1= find(strcmp(namecp,'C1'));
                    model.gammak(indH2)=-1.0;
                    model.gammak(indH2O)=0.5;
                    model.gammak(indCO2)=-0.25;
                    model.gammak(indC1)=0.25;
                end
              %==================SDS MODIF=========================
                 
            end
            model.compFluid = compFluid;
                
            eos = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
              
            model.EOSModel = eos;
            % Set nonlinear tolerance
            model.nonlinearTolerance = 1e-3;
            % Merge in additional optjions
%             model = merge_options(model, varargin{:});
            % Check that we have a valid thermal formulation
            assert(any(strcmpi(model.bacterialFormulation, {'bacterialmodel'})), ...
                'BioChemsitryModel supports currently only one micro-organism');
            % Set output state functions
            model.OutputStateFunctions = {'ComponentTotalMass'};
            if model.bacteriamodel
                model.FlowDiscretization = BiochemicalFlowDiscretization(model);
            end
        end
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@GenericOverallCompositionModel(model);
        end
        
        function model = setupOperators(model, G, rock, varargin)
            % Set up operators, potentially accounting for dynamic
            % transmissibilites

            % Set rock and grid from model if not provided
            if nargin < 3, rock = model.rock; end
            if nargin < 2, G = model.G;       end

            drock = rock;
            if model.dynamicFlowTrans()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                drock = rock;
                nbact0 = 50*10^7;
                drock.perm = rock.perm(1*barsa, nbact0);
            end

            if model.dynamicFlowPv()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                drock = rock;
                nbact0 = 50*10^7;
                drock.poro = rock.poro(0.2, nbact0);
            end
            % Let reservoir model set up operators
            model = setupOperators@ReservoirModel(model, G, drock, varargin{:});
            model.rock = drock;
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate model to see if it is ready for simulation
            % Check that we have a facility model
            if model.bacteriamodel
                if isempty(model.FacilityModel) ...
                        || ~isa(model.FacilityModel, 'BiochemistryGenericFacilityModel')
                    model.FacilityModel = BiochemistryGenericFacilityModel(model);
                end
            else
                if isempty(model.FacilityModel) ...
                        || ~isa(model.FacilityModel, 'GenericFacilityModel')
                    model.FacilityModel = GenericFacilityModel(model);
                end
            end
            model = validateModel@GenericOverallCompositionModel(model, varargin{:});
        end


        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
            % Set up state function groupings for geothermal simulation


            %model.PVTPropertyFunctions = []; % make sure this ir reset
            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});
            fluxprops = model.FlowDiscretization;
            pvtprops  = model.PVTPropertyFunctions;
            flowprops = model.FlowPropertyFunctions;
            if model.bacteriamodel
                flowprops = flowprops.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                flowprops = flowprops.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                flowprops = flowprops.setStateFunction('BactConvRate', BactConvertionRate(model));
            end
            model.PVTPropertyFunctions  = pvtprops;
            model.FlowPropertyFunctions = flowprops;
            model.FlowDiscretization    = fluxprops;

        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Validate state and check if it is ready for simulation

            % Let parent model do it's thing
            state = validateState@GenericOverallCompositionModel(model, state);
            % Set component overall mass fractions if they are not given
            ncomp = length(model.EOSModel.getComponentNames());
            if isfield(state, 'components')
                assert(all(size(state.components) == [model.G.cells.num, ncomp]));
            else
                state.components = zeros(model.G.cells.num, ncomp);
                state.components(:,1) = 1;
            end

            if model.bacteriamodel

                if ~isfield(state, 'nbact')
                    % Set temperature if it is not given
                    nbact0 = 10^6;
                    state.nbact = repmat(nbact0, model.G.cells.num, 1);
                end

            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            [vars, names, origin] = getPrimaryVariables@GenericOverallCompositionModel(model, state);

            if model.bacteriamodel
                nbact = model.getProps(state, 'bacteriamodel');
                vars = [vars, {nbact}];
                names = [names, {'bacteriamodel'}];
                origin = [origin, {class(model)}];
            end
        end

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
        % Get equations from AD states

            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);

             [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
             comps = cellfun(@(x, y) {x, y}, X(:, model.getLiquidIndex), X(:, model.getVaporIndex), 'UniformOutput', false);
            
            
            
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end

            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                             pressures, sat, mob, rho, ...
                                                             {}, comps, ...
                                                             drivingForces);
            if ~isempty(model.FacilityModel)
                % Add sources
                src = model.FacilityModel.getComponentSources(state);
                eqs = model.insertSources(eqs, src);
            end
            if model.bacteriamodel

                src_rate = model.FacilityModel.getProps(state, 'BactConvRate');
                for i = length(eqs)
                    if ~isempty(src_rate{i})                    
                        eqs{i} = eqs{i} -src_rate{i};
                    end
                end

                [beqs, bflux, bnames, btypes] = model.FlowDiscretization.bacteriaConservationEquation(model, state, state0, dt);             
                fd = model.FlowDiscretization;
                src_growthdecay = model.FacilityModel.getBacteriaSources(fd, state, state0, dt);
                beqs{1} = beqs{1} - src_growthdecay;
                %  treat diffusion separately
                if any(model.bDiffusionEffect > 0)
                    beqs{1} = model.operators.AccDiv(beqs{1}, bflux{1});
                end
            else
                [beqs, bnames, btypes] = deal([]);
            end
            % Concatenate
            eqs   = [eqs, beqs];
            names = [names, bnames];
            types = [types, btypes];

            if ~isempty(model.FacilityModel)
                % Get facility equations
                [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
                % Concatenate
                eqs   = [eqs  , weqs  ];
                names = [names, wnames];
                types = [types, wtypes];
            end
            % Add in boundary conditions
%             eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, drivingForces);
            
        end
       
        function forces = validateDrivingForces(model, forces, varargin)
            forces = validateDrivingForces@OverallCompositionCompositionalModel(model, forces, varargin{:});
            forces = validateCompositionalForces(model, forces, varargin{:});
        end 

%         function state = initStateAD(model, state, vars, names, origin)        
%             state = initStateAD@GenericOverallCompositionModel(model, state, vars, names, origin);
%         end
        %-----------------------------------------------------------------%
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
        % Get values for convergence check
        
            [v_eqs, tolerances, names] = getConvergenceValues@ReservoirModel(model, problem, varargin{:});            
        end
        
       
        %-----------------------------------------------------------------%
 function [fn, index] = getVariableField(model, name, varargin)
        switch(lower(name))
          case {'nbact', 'bacteriamodel'} %Bacteria model
                index = ':';
                fn = 'nbact';
          %========================================
            otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@OverallCompositionCompositionalModel(model, name, varargin{:});
        end
    end
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
        % Get names of the fluid components
        
            names  = getComponentNames@ReservoirModel(model);
            names  = horzcat(names,  model.compFluid.names);
            
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dz, drivingForces)
            % Update the state based on increments of the primary values

%             switch model.bacterialFormulation
%                 case 'bacterialmodel'
%                     % Temperature
%                     state = model.updateStateFromIncrement(state, dz, problem, 'bacterialmodel', model.dnMaxRel, model.dnMaxAbs);
%                     state = model.capProperty(state, 'bacterialmodel', model.minimumTemperature, model.maximumTemperature);
%             end
% 
% %             [state, report] = updateState@ReservoirModel(model, state, problem, dz, drivingForces);
                        
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dz, drivingForces);


        end
        
        
        %-----------------------------------------------------------------%
        function isDynamic = dynamicFlowTrans(model)
        % Get boolean indicating if the fluid flow transmissibility is
        % dynamically calculated
        
           isDynamic = isa(model.rock.perm, 'function_handle');
           
        end
        
        %-----------------------------------------------------------------%
        function isDynamic = dynamicFlowPv(model)
        % Get boolean indicating if the fluid flow porevolume is
        % dynamically calculated
        
           isDynamic = isa(model.rock.poro, 'function_handle');
           
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
