classdef ComponentDispersionPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties
        
    end

    methods
        function gp = ComponentDispersionPhaseFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn('PhaseFlux', 'FlowDiscretization'); %flux linked to Darcy velocity
            gp = gp.dependsOn('pressure', 'state');
            gp = gp.dependsOn('mobility', 'state');
           
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
               [Darcy_flux] = model.getProp(state, 'PhaseFlux'); %PhaseFlux linked to Darcy velocity
               op = model.operators;
               [pf]=model.getProp(state, 'pressure');
              
               avg = model.operators.faceAvg;
               Face_poro= avg(model.rock.poro); %average porosity on faces

               interior_faces = find(all(model.G.faces.neighbors ~= 0, 2));
               interior_areas = model.G.faces.areas(interior_faces);
               

               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
               %Grad_face_p=op.Grad(pf);
               %Grad_cell_p=vectorCellGradient(model,Grad_face_p);
 
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
                       u_ph=Darcy_flux{ph}./(interior_areas.*Face_poro);
                       rho_ph=op.faceUpstr(Darcy_flux{ph}, rho{ph});
                       %%fprintf('Darcy_flux{ph}: %16.8f ,  %16.8f \n',min(Darcy_flux{ph}.val),max(Darcy_flux{ph}.val));
                       %%fprintf('u_ph: %16.8f ,  %16.8f \n',min(u_ph.val),max(u_ph.val));
                           
                       if (ph==L_ix) 
                           D_disp = rho_ph.*model.alphaw_long.*(abs(u_ph)+1.e-12);
                          % fprintf('u_L: %16.8f ,  %16.8f \n',min(u_ph.val),max(u_ph.val));
                           %fprintf('D_disp L: %16.8f ,  %16.8f \n',min(D_disp.val),max(D_disp.val));
                           J{c, ph} = - D_disp.*op.Grad(xc);
                           
                       elseif (ph==V_ix) 
                            %fprintf('uG: %8.4f ,  %8.4f \n',min(abs(uG).val),max(abs(uG).val));
                           % fprintf('u_G: %16.8f ,  %16.8f \n',min(u_ph.val),max(u_ph.val));
                           D_disp = rho_ph.*model.alphag_long.*(abs(u_ph)+1.e-12);
                           %fprintf('D_disp G: %16.8f ,  %16.8f \n',min(D_disp.val),max(D_disp.val));
                           J{c, ph} = -D_disp.*op.Grad(yc);
                           
                       end
                   end
               end
            end
            
        end
    end
end
                  