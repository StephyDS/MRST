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
         if model.bacteriamodel && model.liquidPhase
             
             ncomp = model.getNumberOfComponents;

             Psigrowth = prop.getEvaluatedDependencies(state, 'PsiGrowthRate'); 
             Y_H2 = model.Y_H2;
             gammak =model.gammak;
             qbiot_temp =  Psigrowth./Y_H2;
             qbiot =cell(ncomp,1);
             
             for c = 1:ncomp            
                qbiot{c} = gammak(c).*qbiot_temp;
            end
         end         
        end
    end
end