function Jinternfaces = projectCellFluxToFacesAD(model, Jcell)
% This function projects a vectorial flux per cell on the internal faces.
%
% SYNOPSIS:
%   Jinternfaces = projectCellFluxToFacesAD(model, Jcell)
%
% PARAMETERS:
%   model - Compositional model, or the reservoir grid itself.
%   Jcell{d} - d composante of the cell flux [Nc×1] (ADI/double)
% RETURNS:
%   Jinternfaces - oriented normal flux on the internal faces (nIF×1).
%
% Author: [Stéphanie Delage Santacreu]
% Date: [16/09/2025]
% Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
% ---------------------------------------------------------------------------

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