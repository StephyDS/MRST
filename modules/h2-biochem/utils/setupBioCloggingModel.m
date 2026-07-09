function [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc, cp, clogModel)
% setupBioCloggingModelMulti -- Add multi-species bio-clogging effects
%
% PARAMETERS:
%   nbact0 - Cell array of initial bacterial concentrations (e.g., {nbact1_0, nbact2_0})
%   nc     - Vector of characteristic concentrations matching the cell array length [nc1, nc2]
%   cp     - Vector of clogging strengths matching the cell array length [cp1, cp2]

if nargin < 5
    clogModel = true;
end

poro0 = model.rock.poro;
perm0 = model.rock.perm(:, 1);

if clogModel
    % 1. Compute the cumulative initial "scale" factor across all species
    num_species = numel(nbact0);
    scale_sum = 0;
    for i = 1:num_species
        scale_sum = scale_sum + cp(i) * (nbact0(:,i) ./ nc(i)).^2;
    end
    scale = 1 + scale_sum;

    % 2. Define the dynamic pore volume multiplier for multiple species
    % pvMult_nbact expects a cell array of concentration vectors
    pvMult_nbact = @(nbact_cell) 1 ./ (1 + scale .* evalCumulativeClog(nbact_cell, nc));

    % 3. Assign porosity handles
    % --- FIX: accept varargin to handle expanded cell array from evaluateFluid ---
    model.fluid.pvMultR = @(p, varargin) pvMult_nbact(varargin);   % <-- only line changed

    poroFun = @(p, nbact) poro0 .* pvMult_nbact(nbact);
    model.rock.poro = poroFun;

    % 4. Define permeability update function (Kozeny–Carman)
    tauFun = @(p, nbact) ((1 - poro0) ./ (1 - poroFun(p, nbact))).^2 .* ...
        (poroFun(p, nbact) ./ poro0).^3;
    permFun = @(p, nbact) perm0 .* tauFun(p, nbact);
    model.rock.perm = permFun;
else
    model.rock.poro = poro0;
    model.rock.perm = perm0;
    model.fluid.pvMultR = @(varargin) 1;
end

end

% Helper function to loop over cell contents during simulation solver steps
function total_clog = evalCumulativeClog(nbact_cell, nc_vec)
total_clog = 0;
for i = 1:numel(nbact_cell)
    total_clog = total_clog + (nbact_cell{i} ./ nc_vec(i)).^2;
end
end
