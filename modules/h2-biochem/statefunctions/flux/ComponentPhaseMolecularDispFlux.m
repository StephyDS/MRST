% =========================================================================
% CLASS: ComponentPhaseMolecularDispFlux
% =========================================================================
% DESCRIPTION:
% This StateFunction calculates the hydrodynamic (mechanical) dispersion
% flux of components within a phase across grid faces. It is designed for
% compositional, biological, or hydrogen-transport models where dispersion
% is a dominant spreading mechanism.
%
% PHYSICAL MODEL: Bear-Scheidegger (1972) / Bear & Bachmat (1990)
% The dispersion tensor D is defined as:
%    D = alpha_T * |v| * I + (alpha_L - alpha_T) * (v * v^T) / |v|
% where 'v' is the interstitial PORE velocity (Darcy velocity / porosity).
%
% The resulting molar/mass flux across a face is computed as:
%    J_{i,alpha} = - rho * phi * A * (n^T * D * n) * Grad(z_{i,alpha})
%
% NOTE ON TPFA GRIDS:
% In standard MRST Two-Point Flux Approximation (TPFA), flow is assumed
% strictly orthogonal to the faces. Consequently, the velocity is parallel
% to the normal vector, causing the Transversal term (alpha_T) to drop out
% of the 1D face projection. The full tensor math is retained here to
% support higher-order flux methods (e.g., MPFA) if used in the future.
% =========================================================================

classdef ComponentPhaseMolecularDispFlux < StateFunction
    properties
        minPorosity = 1e-12;
        % Default values if not provided by model
        % alphaL_water = 5.0e-2; alphaT_water = 5.0e-3;
        % alphaL_gas   = 1.5e-1; alphaT_gas   = 1.5e-2;
        alphaL_water = 1.0e-1; alphaT_water = 1.0e-2;
        alphaL_gas   = 3.0e-1; alphaT_gas   = 3.0e-2;
    end

    methods
        function sf = ComponentPhaseMolecularDispFlux(model, varargin)
            sf@StateFunction(model, varargin{:});

            % Register dependencies for Automatic Differentiation (AD)
            sf = sf.dependsOn('Density', 'PVTPropertyFunctions');
            sf = sf.dependsOn({'s', 'x', 'y', 'pressure', 'mobility'}, 'state');

            % DEPEND ON THE NATIVE AD FLUX to ensure Newton Jacobian is fully coupled
            sf = sf.dependsOn('PhaseFlux');

            % Check for dynamic porosity dependencies (e.g., bio-clogging)
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                % Explicitly split string arguments so it finds 'nbact' in 'state'
                sf = sf.dependsOn('nbact', 'state');
            end
        end

        function J = evaluateOnDomain(sf, model, state)
            % If composition isn't initialized yet, return zero flux
            if ~isfield(state, 'x')
                [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());
                J = repmat({0}, ncomp, nph);
                return;
            end

            % --- Setup Operators and Geometry ---
            op = model.operators;
            G  = model.G;
            dim = G.griddim;
            [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());

            % --- Fetch Properties ---
            rho        = model.getProps(state, 'Density');
            phase_flux = model.getProps(state, 'PhaseFlux'); % Extract true AD flux
            [xc, yc]   = localGetMoleFractions(model, state);

            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                phi = model.rock.poro(p, nbact);
            else
                phi = model.rock.poro;
            end
            phi_f = op.faceAvg(max(phi, sf.minPorosity));

            % --- Geometry Logic (Internal Faces Only) ---
            intFaces = find(op.internalConn);
            N = op.N; ci = N(:,1); cj = N(:,2);
            ccent = G.cells.centroids;
            A_int = G.faces.areas(intFaces);

            % Compute Unit Normals and Distances
            rij = ccent(cj,:) - ccent(ci,:);
            n_vec = G.faces.normals(intFaces, :);
            sgn = sign(sum(n_vec .* rij, 2));
            n_unit = (n_vec ./ A_int) .* sgn;
            dist_ij = max(sum(n_unit .* rij, 2), 1e-12);

            % --- Dispersivity Map ---
            L_ix = model.getLiquidIndex(); V_ix = model.getVaporIndex();
            aL = zeros(1, nph); aT = zeros(1, nph);
            aL(L_ix) = sf.alphaL_water; aT(L_ix) = sf.alphaT_water;
            if V_ix ~= L_ix
                aL(V_ix) = sf.alphaL_gas; aT(V_ix) = sf.alphaT_gas;
            end

            J = cell(ncomp, nph); [J{:}] = deal(0);

            % --- PHASE LOOP (Velocity & Tensor Physics) ---
            for ph = 1:nph
                rho_f = op.faceAvg(ifcell(rho, ph));

                % 1. Extract Native AD Darcy Flux for the phase
                u_vol_full = ifcell(phase_flux, ph);
                % Safeguard: If PhaseFlux returns all faces, filter to internal faces.
                % If it already returns only internal connections, use directly.
                if size(u_vol_full, 1) == G.faces.num
                    u_vol = u_vol_full(intFaces);
                else
                    u_vol = u_vol_full;
                end

                % 2. Darcy Velocity Vector and Norm
                unorm_sq = 1e-6; % smoothing to prevent div by zero in AD
                u_face = cell(1, dim);
                for d = 1:dim
                    u_face{d} = (u_vol ./ A_int) .* n_unit(:, d);
                    unorm_sq  = unorm_sq + u_face{d}.^2;
                end
                unorm = unorm_sq.^0.5;

                % Convert Darcy velocity to Pore velocity for dispersion tensor
                v_pore = unorm ./ phi_f;

                % 3. Tensor Projection (n' * D * n)
                u_dot_n = 0;
                for d = 1:dim
                    u_dot_n = u_dot_n + (u_face{d} ./ unorm) .* n_unit(:, d);
                end

                % n'Dn = alphaT*|v| + (alphaL - alphaT)*|v|*(u_hat . n)^2
                n_dot_Dn = (aT(ph) .* v_pore) + (aL(ph) - aT(ph)) .* v_pore .* (u_dot_n.^2);

                % --- COMPONENT LOOP (Just Gradients) ---
                for c = 1:ncomp
                    z = pick_mole_frac(ph, L_ix, V_ix, xc{c}, yc{c});
                    grad_z_n = op.Grad(z) ./ dist_ij;

                    % Dispersive flux: J = - rho * phi * A * D_n * dz/dn
                    J{c, ph} = - rho_f .* phi_f .* A_int .* n_dot_Dn .* grad_z_n;
                end
            end
        end
    end
end

% --- Helper Functions ---
function val = ifcell(field, ph)
if iscell(field)
    val = field{ph};
else
    val = field(:, ph);
end
end

function z = pick_mole_frac(ph, L_ix, V_ix, x, y)
if ph == L_ix
    z = x;
elseif ph == V_ix
    z = y;
else
    error('Phase mismatch');
end
end

function [xc, yc] = localGetMoleFractions(model, state)
[x, y] = model.getProps(state, 'x', 'y');
if ~iscell(x)
    xc = mat2cell(x, size(x,1), ones(1, size(x,2)));
else
    xc = x;
end

if ~iscell(y)
    yc = mat2cell(y, size(y,1), ones(1, size(y,2)));
else
    yc = y;
end
end