function dispInfo(GV)
   GV2 = GV.surfGrid;
   fprintf('    -------------------------------------------------------------------------\n');
    fprintf(['  | * Info: The dimensions of Cartesian region of GV ',...
        'is ny = %2.0f, nz = %2.0f.    |\n'],  GV2.cartDims(2), GV.layers.num);
    fprintf(['  |         Users should specify the logical indices of ', ...
        'well region in       |\n']);
    fprintf(['  |         y-z plane within the range of   1 < Ind(1) < Ind(2) < %2.0f  ', ...
        '       |\n'], GV2.cartDims(2));
    fprintf(['  |                                         1 < Ind(3) < Ind(4) < %2.0f  ', ...
        '       |\n'], GV.layers.num);
    fprintf('    -------------------------------------------------------------------------\n');
end