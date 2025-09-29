function gradCell = vectorCellGradient(model, faceGrad)
    % Inputs:
    %   model    : MRST model
    %   faceGrad : [nFaces x 1] GenericAD or double, e.g. op.Grad(p)
    %
    % Output:
    %   gradCell : [nCells x dim] vector gradient of scalar field

    %Mesh informations
    G = model.G;
    ndim = G.griddim;
    nCells = G.cells.num;

    %Geometric informations
    cellVolume=G.cells.volumes; % [nCells,1]
    faceNormals = G.faces.normals;      % [nFaces, ndim]
    faceAreas = sqrt(sum(faceNormals.^2, 2));  % [nFaces, 1]
    unitNormals = faceNormals ./ faceAreas;    % [nFaces, ndim]

    %Face-Cell connectivities
    facePos=G.cells.facePos;
    faces = G.cells.faces;


    % Init output
    gradCell = cell(1, ndim);
    for d = 1:ndim
        if isa(faceGrad, 'ADI') || isa(faceGrad, 'GenericAD')
            gradCell{d} = model.AutoDiffBackend.initVariablesAD(zeros(nCells, 1));
        else
            gradCell{d} = zeros(nCells, 1);
        end
    end


    for i = 1:nCells
       facesLoc=faces(facePos(i):facePos(i+1)-1,1); %global numerotation of faces
       nfacesloc=numel(facesLoc);

       %compute the sign of each face: +1 if normal points outward, -1 if inward
       signs = zeros(nfacesloc, 1);
       for j= 1:nfacesloc
           f=facesLoc(j);
           N=model.operators.N(f,:); % cells neighbors of face f
           if N(1)==i
               signs(j)=1;
           elseif N(2)==i
               signs(j)=-1;
           else
               error('Cell %d not connected to face %d.', i, f);
           end
       end

       %compute the contributions for the cellgradient
       for j= 1:nfacesloc
           f=facesLoc(j);
           nf=unitNormals(f,:);
           Af=faceAreas(f);
           for d = 1:ndim
                contrib = signs(j) .* Af .* faceGrad(f) .* nf(d);
                gradCell{d}(i) = gradCell{d}(i) + contrib;
           end
       end
       for d = 1:ndim
            gradCell{d}(i) = gradCell{d}(i) ./ cellVolume(i);
       end
      
    end

end
