classdef ComponentPhaseHydrodynamicDispFlux < StateFunction
    % Unified mechanical dispersion + molecular diffusion flux for each component/phase.
    %
    % J_{i,alpha} = - rho_alpha * phi_f * A_f * D_total * (dz/dn)
    %
    % D_total = D_disp + D_diff, projected onto the face normal.
    %
    % The class checks the model flags
    %   model.molecularDispersion  (true/false)
    %   model.molecularDiffusion   (true/false)
    % and activates the corresponding physics automatically.
    properties
        minPorosity = 1e-12;
        % --- Mechanical dispersion parameters ---
        alphaL_water = 5.0e-2;   % [m]
        alphaT_water = 5.0e-3;
        alphaL_gas   = 1.5e-1;
        alphaT_gas   = 1.5e-2;
    end
    methods
        function sf = ComponentPhaseHydrodynamicDispFlux(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf = merge_options(sf, varargin{:});
            % Dependencies always needed
            sf = sf.dependsOn('Density', 'PVTPropertyFunctions');
            sf = sf.dependsOn({'s', 'x', 'y', 'pressure'}, 'state');
            % Additional dependencies only when the corresponding flag is on
            if isprop(model, 'molecularDispersion') && model.molecularDispersion
                sf = sf.dependsOn('PhaseFlux');
            end
            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                sf = sf.dependsOn({'MolecularTransmissibility'});
            end
            % Dynamic porosity (e.g., bio‑clogging)
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                sf = sf.dependsOn('nbact', 'state');
            end
            sf.label = 'J_{i,\\alpha}^{disp+diff}';
        end
        function J = evaluateOnDomain(sf, model, state)
            if ~isfield(state, 'x')
                [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());
                J = repmat({0}, ncomp, nph);
                return;
            end
            op = model.operators;
            G  = model.G;
            dim = G.griddim;
            [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            nm = model.getPhaseNames();   % for saturation look‑up
            % --- Geometry (internal faces only) ---------------------------------
            intFaces = find(op.internalConn);
            N = op.N;  ci = N(:,1);  cj = N(:,2);
            ccent = G.cells.centroids;
            A_int = G.faces.areas(intFaces);
            n_vec = G.faces.normals(intFaces, :);
            rij = ccent(cj,:) - ccent(ci,:);
            sgn = sign(sum(n_vec .* rij, 2));
            n_unit = (n_vec ./ A_int) .* sgn;
            dist_ij = max(sum(n_unit .* rij, 2), 1e-12);
            % --- Porosity at faces ---------------------------------------------
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                [p, nbact] = model.getProps(state, 'pressure', 'nbact');
                phi = model.rock.poro(p, nbact);
            else
                phi = model.rock.poro;
            end
            phi_f = op.faceAvg(max(phi, sf.minPorosity));

            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                % Retrieve microbial transmissibility
                Tmoldiff = sf.getEvaluatedDependencies(state, 'MolecularTransmissibility');
            end

            % --- Density -------------------------------------------------------
            rho = model.getProps(state, 'Density');
            [xc, yc] = localGetMoleFractions(model, state);
            % --- Dispersivities (if dispersion enabled) ------------------------
            if isprop(model, 'molecularDispersion') && model.molecularDispersion
                phase_flux = model.getProps(state, 'PhaseFlux');
                aL = zeros(1, nph); aT = zeros(1, nph);
                aL(L_ix) = sf.alphaL_water; aT(L_ix) = sf.alphaT_water;
                if V_ix ~= L_ix
                    aL(V_ix) = sf.alphaL_gas; aT(V_ix) = sf.alphaT_gas;
                end
            end
           
            % --- Output initialisation ------------------------------------------
            J = cell(ncomp, nph); [J{:}] = deal(0);
            
            % --- Phase loop ----------------------------------------------------
            for ph = 1:nph
                rho_f = op.faceAvg(ifcell(rho, ph));
                s = model.getProp(state, ['s', nm(ph)]);
                s_f = op.faceAvg(s);
                % >>> DISPERSION PART <<<
                D_disp = 0;   % face‑wise dispersion coefficient
                if isprop(model, 'molecularDispersion') && model.molecularDispersion
                    u_vol_full = ifcell(phase_flux, ph);
                    if size(u_vol_full, 1) == G.faces.num
                        u_vol = u_vol_full(intFaces);
                    else
                        u_vol = u_vol_full;
                    end
                    unorm_sq = 0;   % correct initialization
                    u_face = cell(1, dim);
                    for d = 1:dim
                        u_face{d} = (u_vol ./ A_int) .* n_unit(:, d);
                        unorm_sq  = unorm_sq + u_face{d}.^2;
                    end
                    unorm = (unorm_sq + 1e-12).^0.5;   % smooth sqrt
                    v_pore = unorm ./ (phi_f.*s_f+sf.minPorosity);            % interstitial velocity
                    u_dot_n = 0;
                    for d = 1:dim
                        u_dot_n = u_dot_n + (u_face{d} ./ unorm) .* n_unit(:, d);
                    end
                    % n^T D_disp n
                    D_disp = (aT(ph) .* v_pore) + (aL(ph) - aT(ph)) .* v_pore .* (u_dot_n.^2);
                end
                % >>> DIFFUSION PART <<<   
                D_diff=0;
                % --- Component loop --------------------------------------------
                for c = 1:ncomp
                    z = pick_mole_frac(ph, L_ix, V_ix, xc{c}, yc{c});
                    grad_z_n = op.Grad(z) ./ dist_ij;   % Δz/Δx
                    if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                       D_diff = Tmoldiff{c,ph};
                    end
                    
                    J{c, ph}=- rho_f .* phi_f .*s_f .* A_int .* D_disp .* grad_z_n ...
                        - rho_f .* D_diff .* op.Grad(z);
                end
            end
        end
    end
end
% ------------------------------------------------------------------------
% Helper functions
% ------------------------------------------------------------------------
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
    if ~iscell(x), xc = mat2cell(x, size(x,1), ones(1, size(x,2))); else xc = x; end
    if ~iscell(y), yc = mat2cell(y, size(y,1), ones(1, size(y,2))); else yc = y; end
end
