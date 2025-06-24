classdef ComponentDispersionPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties
        
    end

    methods
        function gp = ComponentDispersionPhaseFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn('pressure', 'state');
            gp = gp.dependsOn('FaceMobility', 'FlowDiscretization');
            gp = gp.dependsOn('x', 'state');
            gp = gp.dependsOn('y', 'state');
             
            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)            
            ncomp = model.getNumberOfComponents;            
            nph = model.getNumberOfPhases;
            J = cell(ncomp, nph);
            J = cellfun(@(x) 0, J, 'UniformOutput', false);      
            if isfield(state,'x')
               rho = prop.getEvaluatedExternals(model, state, 'Density'); 
               [p, mob_face] = model.getProps(state, 'pressure', 'FaceMobility'); %SDS modeif
               op = model.operators;
               Gradp=op.Grad(p);
   
               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
               
 
               for c = 1:ncomp
                   if iscell(state.x)
                       xc = state.x{c};
                   else
                       xc = state.x(c);
                   end
                   if iscell(state.y)
                       yc = state.y{c};
                   else
                       yc = state.y(c);
                   end
                   for ph = 1:nph
                       u_ph=-op.T.*mob_face{ph}.*Gradp;
                       rho_ph=op.faceUpstr(u_ph, rho{ph});
                       
                       if (ph==L_ix) %&& isfield(state, 'x') && numel(state.x) >= c
                           D_disp = rho_ph.*model.alphaw_long.*(abs(u_ph)+1.e-12);
                           %fprintf('uL: %8.4f ,  %8.4f \n',min(abs(uL).val),max(abs(uL).val));
                           J{c, ph} = - D_disp.*op.Grad(xc);
                           %fprintf('Liq J{c, ph}: %8.4f  \n',J{c, ph});
                           
                       elseif (ph==V_ix) %&& isfield(state, 'y') && numel(state.y) >= c
                            %fprintf('uG: %8.4f ,  %8.4f \n',min(abs(uG).val),max(abs(uG).val));
                           D_disp = rho_ph.*model.alphag_long.*(abs(u_ph)+1.e-12);
                           J{c, ph} = -D_disp.*op.Grad(yc);
                           %fprintf('Gas J{c, ph}: %8.4f  \n',J{c, ph});
                           
                       end
                   end
               end
            end
            
        end
    end
end
                  