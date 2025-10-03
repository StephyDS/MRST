function gradCell = vectorCellGradient(model, faceGrad)
%VECTORCELLGRADIENT  ∇p par cellule à partir de op.Grad(p) sur faces
% - Compatible ADI/GenericAD
% - Gère uniquement les faces internes (cohérent avec op.Grad)
% - Retourne un cell-array {1×dim}, 1 colonne par composante du gradient

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




% function gradCell = vectorCellGradient(model, faceGrad)
%     G   = model.G;
%     op  = model.operators;
%     ndim = G.griddim;
% 
%     faceNormals = G.faces.normals;          % [nFaces × ndim]
%     faceAreas   = G.faces.areas;            % [nFaces, 1]
%     unitNormals = faceNormals ./ faceAreas; % [nFaces × ndim]
% 
%     % Initialiser résultat comme cellule
%     gradCell = cell(1, ndim);
% 
%     for d = 1:ndim
%         % flux artificiel = faceGrad .* (aire * normale_d)
%         flux_d = faceGrad .* (faceAreas .* unitNormals(:,d));
% 
%         % divergence par cellule = somme des flux entrants/sortants
%         div_d  = op.Div(flux_d);
% 
%         % gradient = divergence / volume cellule
%         gradCell{d} = div_d ./ G.cells.volumes;
%     end
% end
% 
% 
% 
% 

% function gradCell = vectorCellGradient(model, faceGrad)
%     % Inputs:
%     %   model    : MRST model
%     %   faceGrad : [nFaces x 1] GenericAD or double, e.g. op.Grad(p)
%     %
%     % Output:
%     %   gradCell : [nCells x dim] vector gradient of scalar field
% 
%     %Mesh informations
%     G = model.G;
%     ndim = G.griddim;
%     nCells = G.cells.num;
% 
%     %Geometric informations
%     cellVolume=G.cells.volumes; % [nCells,1]
%     faceNormals = G.faces.normals;      % [nFaces, ndim]
%     faceAreas = sqrt(sum(faceNormals.^2, 2));  % [nFaces, 1]
%     unitNormals = faceNormals ./ faceAreas;    % [nFaces, ndim]
% 
%     %Face-Cell connectivities
%     facePos=G.cells.facePos;
%     faces = G.cells.faces; %toutes les faces (interieur+frontiere)
% 
% 
%     % Init output
%     gradCell = cell(1, ndim);
%     for d = 1:ndim
%         if isa(faceGrad, 'ADI') || isa(faceGrad, 'GenericAD')
%             gradCell{d} = model.AutoDiffBackend.initVariablesAD(zeros(nCells, 1));
%         else
%             gradCell{d} = zeros(nCells, 1);
%         end
%     end
% 
% 
%     for i = 1:nCells
%        facesLoc=faces(facePos(i):facePos(i+1)-1,1); %global numerotation of faces
%        nfacesloc=numel(facesLoc);
% 
%        %compute the sign of each face: +1 if normal points outward, -1 if inward
%        signs = zeros(nfacesloc, 1);
%        for j= 1:nfacesloc
%            f=facesLoc(j);
%            N=G.faces.neighbors(f,:); % cells neighbors of face f
%            if N(1)==i
%                signs(j)=1;
%            elseif N(2)==i
%                signs(j)=-1;
%            else
%                error('Cell %d not connected to face %d.', i, f);
%            end
%        end
% 
%        %compute the contributions for the cellgradient
%        for j= 1:nfacesloc
%            f=facesLoc(j);
%            nf=unitNormals(f,:);
%            Af=faceAreas(f);
%            for d = 1:ndim
%                contribVec = gradCell{d}*0;
%                contribVec(i) = signs(j) .* Af .* faceGrad(f) .* nf(d);
%                gradCell{d} = gradCell{d} + contribVec;
%                gradCell{d}(i) = gradCell{d}(i) + contrib;
%            end
%        end
%        for d = 1:ndim
%             gradCell{d}(i) = gradCell{d}(i) ./ cellVolume(i);
%        end
% 
%     end
% 
% end
