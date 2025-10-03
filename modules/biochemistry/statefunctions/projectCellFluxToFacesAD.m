function Jinternfaces = projectCellFluxToFacesAD(model, Jcell)
% PROJECTCELLFLUXTOFACES_AD
% Projette un flux vectoriel par cellule J (cell-array {1×dim})
% sur les faces internes -> flux normal orienté (nIF×1), AD-safe.
%
% Entrée:
%   Jcell{d}  : composante d du flux par cellule [Nc×1] (ADI/double)
% Sortie:
%   Fint      : flux normal orienté par face interne [nIF×1] (ADI/double)

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
        tmp = op.faceAvg(Jcell{d});                % tailles possibles: nF×1 ou nIF×1
        if size(tmp.val,1) ~= size(N,1)   % vérifier longueur réelle du vecteur
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