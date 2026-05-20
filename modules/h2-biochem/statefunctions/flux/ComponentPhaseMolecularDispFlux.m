classdef ComponentPhaseMolecularDispFlux < StateFunction
    % Molecular dispersion flux 
    % Jdisp_{i,alpha} = - Kdisp_{i,alpha} * Grad(z_{i,alpha})
    % Kdisp_{i,alpha} = alpha_T * norm(v_Darcy) + 
    % (alpha_L-alpha_T)*v_Darcy * Transpose(v_Darcy)/norm(v_Darcy)
    %
    % alphaL is the Longitudinal dispersivity coefficient
    % alphaT is the Transversal dispersivity coefficient
    % v_Darcy is the Darcy velocity
    %
    % Sources:
    % 1/Impact of hydrogen on the hydrodynamic and bio-chemical
    % behavior, Hagemann & al., DOI 10.1007/s10596-015-9515-6
    % 2/A. E. Scheidegger, 
    % The physics of ﬂow through porous media, The MacMillan Company, 
    % New York, 1957.
    %3/ Bear-Scheiddeger model, and J. Bear and Y. Bachmat, Introduction to modeling 
    % of transport phenomena in porous media, Springer-Verlag,
    % New York, 1990.
   properties
         minPorosity = 1e-12;      % Minimum porosity for numerical stability
     
    end


    methods
        function sf = ComponentPhaseMolecularDispFlux(model, varargin)
            sf@StateFunction(model, varargin{:});

            % Core dependencies
            sf = sf.dependsOn('Density', 'PVTPropertyFunctions');

            % Porosity in h2-biochem is a FlowDiscretization statefunction (BactPorosity)
            % Do not depend on Porosity from FlowDiscretization. Porosity may be static
            % (numeric rock.poro) or dynamic (function handle rock.poro(p, nbact)).

            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                sf = sf.dependsOn('nbact', 'state');
            end

            % State deps
            sf = sf.dependsOn('s', 'state');
            sf = sf.dependsOn('x', 'state');
            sf = sf.dependsOn('y', 'state');
            sf = sf.dependsOn('pressure', 'state');
            sf = sf.dependsOn('mobility', 'state');

            sf.label = 'J_{i,\\alpha}';
        end

        function J = evaluateOnDomain(sf, model, state)
            if ~isfield(state, 'x')
                ncomp = model.getNumberOfComponents();
                nph   = model.getNumberOfPhases();
                J = cell(ncomp, nph);
                [J{:}] = deal(0);
                return;
            end
            J = localMolecularDispersionFluxes(sf, model, state);
        end
    end
end

function J = localMolecularDispersionFluxes(sf, model, state)
op    = model.operators;
ncomp = model.getNumberOfComponents();
nph   = model.getNumberOfPhases();

% --- Porosity (dimensionless) from rock (static or dynamic)

if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
    [p, nbact] = model.getProps(state, 'pressure', 'nbact');
    phi = model.rock.poro(p, nbact);
else
    phi = model.rock.poro;
end
phi_safe = max(phi, sf.minPorosity);

% --- Density
rho = sf.getEvaluatedExternals(model, state, 'Density');

% Mole fractions
[xc, yc] = localGetMoleFractions(model, state);

% Initialize
J = cell(ncomp, nph);
[J{:}] = deal(0);
L_ix = model.getLiquidIndex();
V_ix = model.getVaporIndex();

dim=model.G.griddim;
ncells=model.G.cells.num;
alphaL=zeros(2,1);%Longitudinal dispersivity coefficient
alphaT=zeros(2,1);%Transversal dispersivity coefficient
alphaL(L_ix)=5.e-2; %1.e-2; %5.e-3;%Longitudinal dispersivity coefficient in water phase,  m (0.01->1m)
alphaT(L_ix)=5.e-3;%1.e-3;%5.e-4;%Transversal dispersivity coefficient in water phase,  m
alphaL(V_ix)=1.5e-1; %1.e-2;  %Longitudinal dispersivity coefficient in gas phase,  m (0.1->5m)
alphaT(V_ix)=1.5e-2;%1.e-3;%Transversal dispersivity coefficient in gas phase,  m
% Gas properties
[pf, mob] = model.getProps(state, 'pressure', 'mobility');
Kcell= model.rock.perm; %permeabilities in cells
Grad_cell_p=vectorCellGradient(model,op.Grad(pf));
uDarcy_ph = cell(1, dim);
for d = 1:dim
    uDarcy_ph{d} = Grad_cell_p{d}*0; % init AD-safe
end

for c = 1:ncomp
    for ph = 1:nph
        if iscell(rho)
            rho_ph = rho{ph};
        else
            rho_ph = rho(:,ph);
        end
        if isvector(Kcell) && (numel(Kcell) == ncells)
            % (1) Isotrope scalaire : u = -lam*K*grad p
            for d = 1:dim
                %uDarcy_ph{d} = -mob{ph}.* Kcell.*Grad_cell_p{d}; %without gravity
                uDarcy_ph{d} = -mob{ph}.* Kcell.*(Grad_cell_p{d}-rho_ph.*model.gravity(d));
            end

        elseif ismatrix(Kcell) && all(size(Kcell) == [ncells, dim])
            % (2) Anisotrope diagonal : u_d = -lam*K(:,d)*grad_p_d
            for d = 1:dim
                %uDarcy_ph{d} = -mob{ph}.*Kcell(:,d).*Grad_cell_p{d}; %without gravity
                uDarcy_ph{d} = -mob{ph}.*Kcell(:,d).*(Grad_cell_p{d}-rho_ph.*model.gravity(d));
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

        coef1=alphaT(ph).*unorm;
        coef2=alphaL(ph).*unorm-coef1;

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

        Jcell_phc = cell(1, dim);
        for d = 1:dim
            Jcell_phc{d}=uDarcy_ph{d}.*0; % init AD-safe
        end

        if (ph==L_ix)
            Grad_cell_c=vectorCellGradient(model,op.Grad(xc{c}));%calcul gradient par cellule
        elseif (ph==V_ix)
            Grad_cell_c=vectorCellGradient(model,op.Grad(yc{c})); %calcul gradient par cellule
        end
        for di = 1:dim
            for dj=1:dim
                Jcell_phc{di}=Jcell_phc{di}-rho{ph}.*phi_safe.*Ddisp{di,dj}.*Grad_cell_c{dj}; % init AD-safe
            end
        end
        %projection du flux dispersif par cellule(vectoriel)
        %sur les faces(scalaire) suivant la normale sortante
        Jvect_face_phc=projectCellFluxToFacesAD(model, Jcell_phc);
        J{c, ph} =J{c, ph}+Jvect_face_phc;
    end
end
end


function [xc, yc] = localGetMoleFractions(model, state)
if isfield(state, 'rs') || isfield(state, 'rv')
    error('Black-oil diffusion with rs/rv is not supported here.');
end
[x, y] = model.getProps(state, 'x', 'y');
if iscell(x)
    xc = x;
else
    xc = mat2cell(x, size(x, 1), ones(1, size(x, 2)));
end
if iscell(y)
    yc = y;
else
    yc = mat2cell(y, size(y, 1), ones(1, size(y, 2)));
end
end


function gradCell = vectorCellGradient(model, faceGrad)
G  = model.G;
op = model.operators;

% Dimension spatiale
dim = size(G.faces.normals, 2);

% --- Sélection des faces internes pour être compatible avec op.Grad ---
% internalConn est [nF×1] logique ; true pour faces internes
intFaces = find(op.internalConn);
% On aura: length(intFaces) == length(faceGrad) == size(op.N,1)

% --- Normales unitaires alignées avec la convention de N(:,1)->N(:,2) ---
n_all   = G.faces.normals(intFaces, :);    % [nIF × dim]
A_all   = G.faces.areas(intFaces);         % [nIF × 1]
n_unit  = n_all ./ A_all;                  % [nIF × dim]

% Aligne l’orientation des normales avec la direction i->j de N
N   = op.N; % [nIF × 2]
ci  = G.cells.centroids(N(:,1), :);
cj  = G.cells.centroids(N(:,2), :);
rij = cj - ci;                             % vecteur de i vers j
sgn = sign(sum(n_unit .* rij, 2));         % +1 si déjà aligné, -1 sinon
n_unit = n_unit .* sgn;                    % oriente correctement
% (A_all est positif, n_unit est maintenant orienté comme N)

% --- Calcul du gradient par composante via une "divergence" pondérée ---
% Idée: (∇p)_d ≈ (1/Vc) * Σ_f  (A_f * n_f,d * (∂p/∂n)_f)  avec la même
% convention d’orientation que op.Div/op.N.
gradCell = cell(1, dim);
for d = 1:dim
    % flux artificiel pour la composante d : A * n_d * (∂p/∂n)_f
    flux_d = faceGrad .* (A_all .* n_unit(:, d));   % [nIF × 1] AD/double

    % op.Div somme correctement (avec signes) sur les faces internes vers cellules
    div_d  = op.Div(flux_d);                        % [nC × 1] AD/double

    % Normalisation par volume cellule → (∇p)_d
    gradCell{d} = div_d ./ G.cells.volumes;         % [nC × 1]
end
end

function Jinternfaces = projectCellFluxToFacesAD(model, Jcell)
% This function projects a vectorial flux per cell on the internal faces.
G   = model.G;
op  = model.operators;
dim = G.griddim;

% --- Faces internes uniquement (cohérent avec op.N, op.Grad, etc.)
intFaces = find(op.internalConn);             % indices faces internes (longueur nIF)

% --- Normales *aire* restreintes aux faces internes
n_all  = G.faces.normals(intFaces, :);         % [nIF×dim], contient n_f * A_f
A_int  = G.faces.areas(intFaces);              % [nIF×1]
n_unit = n_all ./ A_int;                       % normales unitaires

% --- Aligner l'orientation des normales avec le sens i->j de N
N        = op.N;                               % [nIF×2] indices cellules (i->j)
ci  = G.cells.centroids(N(:,1), :);
cj  = G.cells.centroids(N(:,2), :);
rij = cj - ci;                                 % vecteur de i vers j
sgn = sign(sum(n_unit .* rij, 2));             % +1 si déjà aligné, -1 sinon
n_oriented = n_all .* sgn;                     % contient (n_f * A_f) orienté comme N

% --- Moyenne cell -> face (AD-safe). Selon backends, faceAvg peut rendre nF ou nIF.
Jvect_face = cell(1, dim);
for d = 1:dim
    tmp = op.faceAvg(Jcell{d});
    if isa(tmp, 'GenericAD') || isa(tmp, 'ADI')
        nTmp = size(tmp.val,1);
    else
        nTmp = size(tmp,1);
    end
    if nTmp ~= size(N,1)
        tmp = tmp(intFaces);
    end
    Jvect_face{d} = tmp;                               % [nIF×1] AD/double
end

% --- Projection: flux normal intégré sur la face (orienté i->j)
Jinternfaces = Jvect_face{1} .* n_oriented(:,1);
for d = 2:dim
    Jinternfaces = Jinternfaces + Jvect_face{d} .* n_oriented(:,d);
end
end


