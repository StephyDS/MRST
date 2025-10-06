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
                avg = model.operators.faceAvg;
                interior_faces = find(all(model.G.faces.neighbors ~= 0, 2));
                interior_areas = model.G.faces.areas(interior_faces);
                dim=model.G.griddim;
                ncells=model.G.cells.num;


                [pf, mob] = model.getProps(state, 'pressure','mobility');    % {lambda_w, lambda_g}
                Face_poro= avg(model.rock.poro); %average porosity on faces
                Kcell= model.rock.perm; %permeabilities in cells


                L_ix = model.getLiquidIndex();
                V_ix = model.getVaporIndex();
                Grad_face_p=op.Grad(pf);
                Grad_cell_p=vectorCellGradient(model,Grad_face_p);

                uDarcy_ph = cell(1, dim);
                for d = 1:dim
                    uDarcy_ph{d} = Grad_cell_p{d}*0; % init AD-safe
                end


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
                        if isvector(Kcell) && (numel(Kcell) == ncells)
                            % (1) Isotrope scalaire : u = -lam*K*grad p
                            for d = 1:dim
                                uDarcy_ph{d} = -mob{ph}.* Kcell.*Grad_cell_p{d};
                            end
                        elseif ismatrix(Kcell) && all(size(Kcell) == [ncells, dim])
                            % (2) Anisotrope diagonal : u_d = -lam*K(:,d)*grad_p_d
                            for d = 1:dim
                                uDarcy_ph{d} = -mob{ph}.*Kcell(:,d).*Grad_cell_p{d};
                            end
                        end
                        unorm=uDarcy_ph{1}.^2;
                        for d = 2:dim
                            unorm=unorm+uDarcy_ph{d}.^2;
                        end
                        unorm=unorm.^0.5;
                        uhat  = cell(1, dim);
                        for d = 1:dim
                            uhat{d} = uDarcy_ph{d} ./ (unorm + 1e-12);
                        end
                        coef1=model.alphaT(ph).*unorm;
                        coef2=model.alphaL(ph).*unorm-coef1;

                        Ddisp=cell(dim,dim); %rho*Ddisp par phase dim*dim cell
                        for di=1:dim
                            for dj=1:dim
                                Ddisp{di,dj}=coef1.*0;
                            end
                        end
                        for di=1:dim
                            Ddisp{di,di}= coef1;
                            for dj=1:dim
                                Ddisp{di,dj}=Ddisp{di,dj}+coef2.*uhat{di}.*uhat{dj};
                            end
                        end
                        % if (ph==V_ix)
                        %     for di=1:dim
                        %         for dj=1:dim
                        %             fprintf('ph=%8.1f, i=%8.1f, j=%8.1f, Ddisp : %16.8f , %16.8f \n',...
                        %                 ph,di,dj,min(Ddisp{di,dj}.val),max(Ddisp{di,dj}.val));
                        %         end
                        %     end
                        % end

                        Jcell_phc = cell(1, dim);
                        for d = 1:dim
                            Jcell_phc{d}=uDarcy_ph{d}.*0; % init AD-safe
                        end

                        if (ph==L_ix)
                            Grad_cell_c=vectorCellGradient(model,op.Grad(xc));%calcul gradient par cellule
                        elseif (ph==V_ix)
                            Grad_cell_c=vectorCellGradient(model,op.Grad(yc)); %calcul gradient par cellule
                        end
                        for di = 1:dim
                            for dj=1:dim
                                Jcell_phc{di}=Jcell_phc{di}-rho{ph}.*Ddisp{di,dj}.*Grad_cell_c{dj}; % init AD-safe
                            end
                        end
                        %projection du flux dispersif (vectoriel par
                        %cellule) sur les faces (scalaire) suivant la
                        %normale sortante
                        Jvect_face_phc=projectCellFluxToFacesAD(model, Jcell_phc);
                        %J{c, ph} =Jvect_face_phc;
                        fprintf('ph= %16.2f, Jvect_face_phc: %16.8f ,  %16.8f \n',ph,...
                            min(Jvect_face_phc.val),max(Jvect_face_phc.val));

                        u_ph=Darcy_flux{ph}./(interior_areas.*Face_poro);
                        rho_ph=op.faceUpstr(Darcy_flux{ph}, rho{ph});

                        if (ph==L_ix)
                            D_disp = rho_ph.*model.alphaL(ph).*(abs(u_ph)+1.e-12);
                             %fprintf('u_L: %16.8f ,  %16.8f \n',min(u_ph.val),max(u_ph.val));
                            %fprintf('D_disp L: %16.8f ,  %16.8f \n',min(D_disp.val),max(D_disp.val));
                            J{c, ph} = - D_disp.*op.Grad(xc);
                            %tmp1=J{c, ph};
                            %fprintf('J{c, ph} L: %16.8f ,  %16.8f \n',min(tmp1.val),max(tmp1.val));

                        elseif (ph==V_ix)
                            %fprintf('uG: %8.4f ,  %8.4f \n',min(abs(uG).val),max(abs(uG).val));
                            % fprintf('u_G: %16.8f ,  %16.8f \n',min(u_ph.val),max(u_ph.val));
                            D_disp = rho_ph.*model.alphaL(ph).*(abs(u_ph)+1.e-12);
                            %fprintf('D_disp G: %16.8f ,  %16.8f \n',min(D_disp.val),max(D_disp.val));
                            J{c, ph} = -D_disp.*op.Grad(yc);
                            tmp2=J{c, ph};
                            fprintf('J{c, ph} G: %16.8f ,  %16.8f \n',min(tmp2.val),max(tmp2.val));

                        end
                    end
                end
            end

        end
    end
end
