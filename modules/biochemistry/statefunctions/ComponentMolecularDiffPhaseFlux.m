classdef ComponentMolecularDiffPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties
    end

    methods
        function gp = ComponentMolecularDiffPhaseFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
             gp = gp.dependsOn('s', 'state');
             gp = gp.dependsOn('x', 'state');

            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)            
            ncomp = model.getNumberOfComponents;            
            nph = model.getNumberOfPhases;
            J = cell(ncomp, nph);
            J = cellfun(@(x) 0, J, 'UniformOutput', false);      
            if isfield(state,'x')
               nm = model.getPhaseNames();
               tau = [1,1];
               rho = prop.getEvaluatedExternals(model, state, 'Density');      
               avg = model.operators.faceAvg;
               poro= model.rock.poro;
               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
               %if iscell(rho)
                   %fprintf('nbphase:=%.2f, min(rho{1}.val)=%.2f \n',nph,min(rho{1}.val));
                   %ddd1=rho{1}.*tau(1).*poro;
                   %ddd2=rho{2}.*tau(2).*poro;
                   %fprintf('nbphase:=%.2f, min(ddd1.val)=%.2f, min(ddd2.val)=%.2f \n',nph,min(ddd1.val),min(ddd2.val));
                   
               %else
               %    fprintf('rho not cell \n');
               %end

               % Define diffusion coefficients in m²/s for liquid and gas phases
               % These are example values, please replace them with actual data as needed
               % Format: [liquid_diff gas_diff] for each component
               for c = 1:ncomp
                   for ph = 1:nph
                       s = model.getProp(state, ['s', nm(ph)]);
                       %D_diff = avg(model.mol_diff(c,ph).*tau(ph).*poro);
                       %fprintf('comp:=%.2f, phase=%.2f  min(x(1))=%.2f,min(x(2))=%.2f \n',...
                       %  c,ph,min(state.x{1}.val),min(state.x{1}.val));
                       %ddd1=avg(s.*rho{ph}.*model.mol_diff(c,ph).*tau(ph).*poro);
                       %ddd2=model.operators.Grad(state.x{1});
                       %ddd3=ddd1.*ddd2;
                       %fprintf('nbphase:=%.2f, size(ddd1.val,1)=%.2f size(ddd3.val,1)=%.2f \n',nph,min(ddd1.val),size(ddd3.val,1));
                       %============SDS MODIF TAU MILLINGTON AND QUIRK MODEL======
                       tau_mq=(s.*poro).^(7/3).*poro.^(-2);
                       D_diff = avg(s.*rho{ph}.*model.mol_diff(c,ph).*tau_mq.*poro);%Millington and Quirk model
                       %D_diff = avg(s.*rho{ph}.*model.mol_diff(c,ph).*tau(ph).*poro);%Millington and Quirk model
                       %============SDS MODIF ======
               
                       if (ph==L_ix) 
                           
                           J{c, ph} = - D_diff.*model.operators.Grad(state.x{c});

                       elseif (ph==V_ix)   
                           %fprintf('VAP:J{c, ph}: c=%.2f, ph=%.2f,mol_diff=%.2f, \n',c,ph,model.mol_diff(c,ph));
                          
                           J{c, ph} = - D_diff.*model.operators.Grad(state.y{c});
                       end
                   end
               end
            end
        end
    end
end
