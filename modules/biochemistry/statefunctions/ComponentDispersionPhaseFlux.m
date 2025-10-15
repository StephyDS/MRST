classdef ComponentDispersionPhaseFlux < StateFunction
    % Flux of each component, in each phase
    % Dispersive flux :Bear–Scheidegger (anisotrope)dispersive model
    % Author: [Stéphanie Delage Santacreu]
    % Date: [16/09/2025]
    % Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
    % ---------------------------------------------------------------------------

    properties

    end

    methods
        function gp = ComponentDispersionPhaseFlux(model, varargin)
            gp@StateFunction(model);
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
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
                op = model.operators;
                dim=model.G.griddim;
                ncells=model.G.cells.num;

                [pf, mob] = model.getProps(state, 'pressure','mobility');    % {lambda_w, lambda_g}
                Kcell= model.rock.perm; %permeabilities in cells
                poro= model.rock.poro;

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
                                %uDarcy_ph{d} = -mob{ph}.* Kcell.*Grad_cell_p{d}; %without gravity
                                uDarcy_ph{d} = -mob{ph}.* Kcell.*(Grad_cell_p{d}-rho{ph}.*model.gravity(d));
                            end

                        elseif ismatrix(Kcell) && all(size(Kcell) == [ncells, dim])
                            % (2) Anisotrope diagonal : u_d = -lam*K(:,d)*grad_p_d
                            for d = 1:dim
                                %uDarcy_ph{d} = -mob{ph}.*Kcell(:,d).*Grad_cell_p{d}; %without gravity
                                uDarcy_ph{d} = -mob{ph}.*Kcell(:,d).*(Grad_cell_p{d}-rho{ph}.*model.gravity(d));
                            end
                        end
                        unorm=uDarcy_ph{1}.^2;
                        for d = 2:dim
                            unorm=unorm+uDarcy_ph{d}.^2;
                        end
                        unorm=unorm.^0.5;
                        %fprintf('ph=%8.1f, unorm : %16.8f , %16.8f \n',...
                        %                 ph,min(unorm.val),max(unorm.val));
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
                                Jcell_phc{di}=Jcell_phc{di}-rho{ph}.*poro.*Ddisp{di,dj}.*Grad_cell_c{dj}; % init AD-safe
                            end
                        end
                        %projection du flux dispersif par cellule(vectoriel)
                        %sur les faces(scalaire) suivant la normale sortante
                        Jvect_face_phc=projectCellFluxToFacesAD(model, Jcell_phc);
                        J{c, ph} =Jvect_face_phc;
                    end
                end
            end

        end
    end
end
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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
