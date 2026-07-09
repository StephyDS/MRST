function [nls, lsolve] = setupOptimizedLinearSolver(model, varargin)
% Setup optimized AMGCL CPR linear solver for compositional biochemical models
%
% SYNOPSIS:
%   [nls, lsolve] = setupOptimizedLinearSolver(model)
%   [nls, lsolve] = setupOptimizedLinearSolver(model, 'complexityLevel', 'high')
%
% PARAMETERS:
%   model           - BiochemistryModel instance
%   complexityLevel - 'low' (advection only), 'medium' (with diffusion),
%                     'high' (with both diffusion & dispersion) [default: 'medium']
%   solverTolerance - Tolerance for linear solver [default: 1e-4]
%   maxNonlinIter   - Max nonlinear iterations [default: 15]
%   cprDamp         - CPR damping factor [default: 0.7]
%
% RETURNS:
%   nls    - NonLinearSolver configured with linear solver
%   lsolve - CPRSolverAD linear solver with optimized AMGCL parameters

    opt = struct('complexityLevel', 'medium', ...
                 'solverTolerance', 1e-4, ...
                 'maxNonlinIter', 15, ...
                 'cprDamp', 0.7);
    opt = merge_options(opt, varargin{:});

    % Select base linear solver
    lsolve = selectLinearSolverAD(model);

    % Configure AMGCL CPR parameters based on problem complexity
    amgclSettings = struct();
    amgclSettings.coarsening = 'aggregation';
    amgclSettings.aggregation = 'legacy';
    amgclSettings.maxLvl = 20;
    amgclSettings.verbosity = 0;
    amgclSettings.cpr_damp = opt.cprDamp;

    switch lower(opt.complexityLevel)
        case 'low'
            % Simple advection-only: aggressive coarsening
            amgclSettings.aggEps = 1e-1;
            amgclSettings.aggBlockSize = 2;
            amgclSettings.eps = 1e-2;
            amgclSettings.maxIter = 50;
            amgclSettings.cpr_maxIter = 3;

        case 'medium'
            % Diffusion included: balanced settings (default)
            amgclSettings.aggEps = 1e-2;
            amgclSettings.aggBlockSize = 1;
            amgclSettings.eps = 1e-3;
            amgclSettings.maxIter = 100;
            amgclSettings.cpr_maxIter = 4;

        case 'high'
            % Both diffusion & dispersion: conservative coarsening
            amgclSettings.aggEps = 5e-3;
            amgclSettings.aggBlockSize = 1;
            amgclSettings.eps = 5e-4;
            amgclSettings.maxIter = 150;
            amgclSettings.cpr_maxIter = 5;

        otherwise
            error('Unknown complexityLevel: %s', opt.complexityLevel);
    end

    % Apply AMGCL settings
    if isprop(lsolve, 'amgclSettings')
        lsolve.amgclSettings = amgclSettings;
    end
    lsolve.tolerance = opt.solverTolerance;

    % Configure nonlinear solver
    nls = NonLinearSolver();
    nls.LinearSolver = lsolve;
    nls.maxIterations = opt.maxNonlinIter;
    nls.verbose = false;
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
