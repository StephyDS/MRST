classdef BactConvertionRate <  StateFunction
    % The bacterial growth rate, given per cell
    properties
    end
    
    methods
        function gp = BactConvertionRate(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PsiGrowthRate'}, 'FlowDiscretization');
            gp.label = 'Q_biot';
        end

        function qbiot = evaluateOnDomain(prop, model, state)
            qbiot = 0;
         if model.ReservoirModel.bacteriamodel && model.ReservoirModel.liquidPhase
             
             ncomp = model.getNumberOfComponents;

             Psigrowth = prop.getEvaluatedDependencies(state, 'PsiGrowthRate'); 
             Y_H2 = model.ReservoirModel.Y_H2;
             gammak =model.ReservoirModel.gammak;
             qbiot_temp =  Psigrowth./Y_H2;
             qbiot =cell(ncomp,1);
             
             for c = 1:ncomp            
                qbiot{c} = gammak(c).*qbiot_temp +0;
            end
         end         
        end
    end
end