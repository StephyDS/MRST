classdef ComponentMolecularDiffPhaseFlux < StateFunction
    % Flux of each component, in each phase
    properties
        
    end

    methods
        function gp = ComponentMolecularDiffPhaseFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');%sds modif
             gp = gp.dependsOn('s', 'state');
             gp = gp.dependsOn('x', 'state');
             gp = gp.dependsOn('y', 'state');
             gp = gp.dependsOn('pressure', 'state');
             gp = gp.dependsOn('T', 'state');

            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)            
            ncomp = model.getNumberOfComponents;            
            nph = model.getNumberOfPhases;
            J = cell(ncomp, nph);
            J = cellfun(@(x) 0, J, 'UniformOutput', false);      
            if isfield(state,'x')
               nm = model.getPhaseNames();
               rho = prop.getEvaluatedExternals(model, state, 'Density'); 
               [p, T] = model.getProps(state, 'pressure', 'temperature'); 
               coeff1=T.^(1.5)./(9.869e-6.*p); %1.e5.*p en atm et p en Pa 
               Molmass=1.e3.*model.compFluid.molarMass;
               SigLJ=model.param_LJ(:,1);
               EpsLJ=model.param_LJ(:,2);

               avg = model.operators.faceAvg;
               poro= model.rock.poro;
               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
               
                interior_faces = find(all(model.G.faces.neighbors ~= 0, 2));
                interior_areas = model.G.faces.areas(interior_faces);
                
               % Define diffusion coefficients in m²/s for liquid and gas phases
               % These are example values, please replace them with actual data as needed
               % Format: [liquid_diff gas_diff] for each component
               sqrtMij=zeros(ncomp,ncomp);
               sqrtEpsij=zeros(ncomp,ncomp);
               Sigij2=zeros(ncomp,ncomp);
               for c = 1:ncomp
                   for cj = 1:ncomp
                       sqrtMij(c,cj)=sqrt(2*Molmass(c)*Molmass(cj)/(Molmass(c)+Molmass(cj)));
                       sqrtEpsij(c,cj)=sqrt(EpsLJ(c)*EpsLJ(cj));
                       Sigij2(c,cj)=0.25*(SigLJ(c)+SigLJ(cj))*(SigLJ(c)+SigLJ(cj));
                   end
               end
               Dij=1.e-4.*0.001858./(sqrtMij.*Sigij2); %en m2/s


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
                       s = model.getProp(state, ['s', nm(ph)]);
                       tau_mq=(s.*poro).^(7/3).*poro.^(-2);
                       
                       if (ph==L_ix) 
                           D_diffl = avg(s.*rho{ph}.*model.mol_diff(c,ph).*tau_mq.*poro);%Millington and Quirk model
                           %J{c, ph} = - D_diffl.*model.operators.Grad(xc);
                           J{c, ph} = - D_diffl.*model.operators.Grad(xc).*interior_areas;

                       elseif (ph==V_ix)
                           %calcul des Dij
                           D_diffij_inv =yc.*0;

                           for cj = 1:ncomp
                               if iscell(state.y)
                                   ycj = state.y{cj};
                               else
                                   ycj = state.y(cj);
                               end
                               if (cj~=c)
                                   D_diffij_inv =  D_diffij_inv  +ycj./Dij(c,cj);
                               end
                           end
                           D_diffij=coeff1./D_diffij_inv;
                           D_diff = avg(s.*rho{ph}.*D_diffij.*tau_mq.*poro);%Millington and Quirk model
                           %J{c, ph} = - D_diff.*model.operators.Grad(yc);
                           J{c, ph} = - D_diff.*model.operators.Grad(yc).*interior_areas;
                       end
                   end
               end
            end
        end
    end
end
