classdef BactConvertionRate <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate', 'PsiDecayRate'}, 'FlowDiscretization');
           %gp = gp.dependsOn({'nbact', 's'}, 'state');
            gp = gp.dependsOn({'PoreVolume', 'Density'}, 'PVTPropertyFunctions');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0;
         if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase
             
             ncomp = model.ReservoirModel.EOSModel.getNumberOfComponents;
             namecp = model.ReservoirModel.EOSModel.getComponentNames();%SDS modif
             indH2=find(strcmp(namecp,'H2'));%SDS modif

             Psigrowth = model.getProps(state, 'PsiGrowthRate');
             rho = model.ReservoirModel.PVTPropertyFunctions.get(model.ReservoirModel, state, 'Density');
             L_ix = model.ReservoirModel.getLiquidIndex();

             if ~iscell(rho)
                 rho = {rho};
             end
             
             compn = model.ReservoirModel.EOSModel.getComponentNames();
             idxH2 = find(strcmp(compn, 'H2'));  % Index of CO2

             Y_H2 = model.ReservoirModel.Y_H2;
             gammak = model.ReservoirModel.gammak;
             nbactMax =model.ReservoirModel.nbactMax;
             nbact0 =model.ReservoirModel.nbact0;
             
             qbiot_temp =  (Psigrowth)./Y_H2./abs(gammak(idxH2));
             qbiot = cell(ncomp,1);
             mc = model.ReservoirModel.EOSModel.CompositionalMixture.molarMass;
             for c = 1:ncomp            
                qbiot{c} = gammak(c).*qbiot_temp.*mc(c).*nbactMax./rho{L_ix} + 0;
             end
         end         
        end
    end
end
