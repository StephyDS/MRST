function gradCell = vectorCellGradient(model, faceGrad)
% This function computes a gradient per cell from a gradient per face (internal faces) .
%
% SYNOPSIS:
%   gradCell = vectorCellGradient(model, faceGrad)
%
% PARAMETERS:
%   model - Compositional model, or the reservoir grid itself.
%   faceGrad - gradient per internal face
% RETURNS:
%   gradCell - gradient per cell.
%
% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
% ---------------------------------------------------------------------------

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

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}
