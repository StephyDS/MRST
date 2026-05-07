classdef ComponentPhaseMolecularDiffFlux < StateFunction
    % Molecular diffusion flux for each component/phase.
    %
    % J_{i,alpha} = - K_{i,alpha} * Grad(z_{i,alpha})
    % K_{i,alpha} = (phi*S_alpha*tau_MQ) * rho_alpha * D_{i,alpha}
    %
    % Millington-Quirk:
    %   (phi*S)*tau_MQ = (phi*S)^(7/3) / phi^2

    methods
        function sf = ComponentPhaseMolecularDiffFlux(model, varargin)
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
            sf = sf.dependsOn('temperature', 'state');
            sf = sf.dependsOn('PhaseFlux', 'FlowDiscretization'); %flux linked to Darcy velocity
            
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
            J = localMolecularDiffusionFluxes(sf, model, state);
        end
    end
end

function J = localMolecularDiffusionFluxes(sf, model, state)
op    = model.operators;
ncomp = model.getNumberOfComponents();
nph   = model.getNumberOfPhases();
nm    = model.getPhaseNames();

% --- Porosity (dimensionless) from rock (static or dynamic)

if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
    if model.bacteriamodel
        [p, nbact] = model.getProps(state, 'pressure', 'nbact');
        nbioreact=model.biochemFluid.nbioreact;
        if nbioreact==1
            if iscell(nbact)
                phi = model.rock.poro(p, nbact{1}); % Apply both modifications
            else
                phi = model.rock.poro(p, nbact); % Apply both modifications
            end
        elseif nbioreact==2
            if iscell(nbact)
                phi = model.rock.poro(p, nbact{1}, nbact{2}); % Apply both modifications
            else
                phi = model.rock.poro(p, nbact(:,1), nbact(:,2)); % Apply both modifications
            end
        end
    else
        p = model.getProp(state, 'pressure');
        phi=model.rock.poro(p);
    end
else
    phi = model.rock.poro;
end
phi_safe = max(phi, 1e-12);

% --- Density
rho = sf.getEvaluatedExternals(model, state, 'Density');

% Mole fractions
[xc, yc] = localGetMoleFractions(model, state);

% Initialize
J = cell(ncomp, nph);
[J{:}] = deal(0);
L_ix = model.getLiquidIndex();
V_ix = model.getVaporIndex();

if (model.molecularDiffusion)
    % Diffusion parameters
    [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(model);
    Dij_ref = localBinaryDiffusionReference(model, paramLJ);

    % Gas scaling
    [p, T] = model.getProps(state, 'pressure', 'temperature');
    Tref = 273.15 + 40;
    pref = atm;
    p_safe = max(p, 1e-8*barsa);
    gasScale = (T./Tref).^1.75 .* (pref./p_safe);

    for ph = 1:nph
        s = model.getProp(state, ['s', nm(ph)]);

        % (phi*S)*tau_MQ
        phiS = phi .* s;
        geom = (phiS).^(7/3) .* (phi_safe).^(-2);

        if ph == L_ix
            for c = 1:ncomp
                Kc = geom .* rho{ph} .* Dliq_ref(c);
                Kf = localFaceHarmonicAvg(op, Kc);
                J{c, ph} = -Kf .* op.Grad(xc{c});
            end

        elseif ph == V_ix
            yAll = model.getProps(state, 'y');
            for c = 1:ncomp
                invDi = 0;
                for j = 1:ncomp
                    if j == c, continue; end
                    yj  = localGetComponentVector(yAll, j);
                    Dij = Dij_ref(c, j) .* gasScale;
                    invDi = invDi + yj ./ max(Dij, 1e-30);
                end
                Di_mix = 1 ./ max(invDi, 1e-30);

                Kc = geom .* rho{ph} .* Di_mix;
                Kf = localFaceHarmonicAvg(op, Kc);
                J{c, ph} = -Kf .* op.Grad(yc{c});
            end
        end
    end   
end

if (model.molecularDispersion)
    [Darcy_flux] = model.getProp(state, 'PhaseFlux'); %PhaseFlux linked to Darcy velocity
    avg = model.operators.faceAvg;
    interior_faces = find(all(model.G.faces.neighbors ~= 0, 2));
    interior_areas = model.G.faces.areas(interior_faces);
    Face_poro= avg(phi_safe); %average porosity on faces

    alphaL=zeros(2,1);%Longitudinal dispersivity coefficient
    alphaL(L_ix)=5.e-3;%1.e-2; %Longitudinal dispersivity coefficient in water phase,  m (0.01->1m)
    alphaL(V_ix)=1.e-2;  %5.e-2; %Longitudinal dispersivity coefficient in gas phase,  m (0.1->5m)

    for c = 1:ncomp
        for ph = 1:nph
            u_ph=Darcy_flux{ph}./(interior_areas.*Face_poro);
            if iscell(rho)
                rho_ph = op.faceUpstr(Darcy_flux{ph}, rho{ph});
            else
                rho_ph = op.faceUpstr(Darcy_flux{ph}, rho(:,ph));
            end
            if (ph==L_ix)
                D_disp = rho_ph.*alphaL(ph).*(abs(u_ph)+1.e-12);
               J{c, ph} = J{c, ph}- D_disp.*op.Grad(xc{c});
               
            elseif (ph==V_ix)
                D_disp = rho_ph.*alphaL(ph).*(abs(u_ph)+1.e-12);
                J{c, ph} = J{c, ph}-D_disp.*op.Grad(yc{c});
              
            end
        end
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

function v = localGetComponentVector(z, j)
if iscell(z)
    v = z{j};
else
    v = z(:, j);
end
end

function Kf = localFaceHarmonicAvg(op, Kc)
if isfield(op, 'N')
    N  = op.N;
    i1 = N(:, 1); i2 = N(:, 2);
    K1 = Kc(i1);  K2 = Kc(i2);
    Kf = 2 .* K1 .* K2 ./ max(K1 + K2, 1e-30);
else
    Kf = op.faceAvg(Kc);
end
end

function [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(model)
namecp = model.compFluid.names();
ncomp  = numel(namecp);

Dliq_ref = 1e-9 * ones(ncomp, 1);
paramLJ  = repmat([3.5, 150.0], ncomp, 1);

db.Dliq = struct('h2',6.44e-9,'c1',2.15e-9,'methane',2.15e-9,'co2',2.72e-9, ...
    'h2o',3.29e-9,'water',3.29e-9,'n2',2.86e-9,'c2',1.72e-9, ...
    'ethane',1.72e-9,'c3',1.43e-9,'propane',1.43e-9,'h2s',2.15e-9, ...
    'nc4',1.15e-9,'butane',1.15e-9);

db.LJ = struct('h2',[2.92,59.7],'c1',[3.758,148.6],'methane',[3.758,148.6], ...
    'co2',[3.996,195.2],'h2o',[2.641,809.1],'water',[2.641,809.1], ...
    'n2',[3.798,71.4],'c2',[4.443,215.7],'ethane',[4.443,215.7], ...
    'c3',[5.118,237.1],'propane',[5.118,237.1],'h2s',[3.60,301.0], ...
    'nc4',[5.206,289.5],'butane',[5.206,289.5]);

for i = 1:ncomp
    key = lower(namecp{i});
    if isfield(db.Dliq, key), Dliq_ref(i) = db.Dliq.(key); end
    if isfield(db.LJ, key),   paramLJ(i,:) = db.LJ.(key);  end
end
end

function Dij_ref = localBinaryDiffusionReference(model, paramLJ)
ncomp = size(paramLJ, 1);
Molmass = 1e3 .* model.compFluid.molarMass; % kg/mol -> g/mol
sigma = paramLJ(:, 1); % Å

Dij_ref = zeros(ncomp);
for i = 1:ncomp
    for j = 1:ncomp
        if i == j
            Dij_ref(i, j) = inf;
            continue;
        end
        sqrtMij   = sqrt(2 * Molmass(i) * Molmass(j) / (Molmass(i) + Molmass(j)));
        sigma_ij2 = 0.25 * (sigma(i) + sigma(j))^2;
        Dij_ref(i, j) = 1e-4 * 0.001858 / (sqrtMij * sigma_ij2);
    end
end
Dij_ref = 0.5*(Dij_ref + Dij_ref.');
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


