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
        % --- Molecular diffusion parameters ---
        Tref = 273.15 + 40;               % [K]
        pref = 101325;                     % [Pa] (1 atm)
        tortuosityExponent = 7/3;          % Millington–Quirk exponent
        gasDiffExponent = 1.5;            % Chapman–Enskog temperature exponent (7/4)
        defaultLiquidDiffusivity = 1e-9;   % [m²/s]
        minDiffusivity = 1e-15;            % reasonable floor for diffusivities
        % Lennard–Jones defaults for unknown components
        defaultSigma   = 3.5;              % [Å]
        defaultEpsilon = 150.0;            % [K]
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
                sf = sf.dependsOn('temperature', 'state');
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
            % --- Diffusion parameters (if enabled) ------------------------------
            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                [p, T] = model.getProps(state, 'pressure', 'temperature');
                p_safe = max(p, 1e-8*barsa);              % ensure pressure > 0
                gasScale = (T./sf.Tref).^sf.gasDiffExponent .* (sf.pref./p_safe);
                [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(sf, model);
                Dij_ref = localBinaryDiffusionReference(sf, model, paramLJ);  % no temperature factor inside
            end
            % --- Output initialisation ------------------------------------------
            J = cell(ncomp, nph); [J{:}] = deal(0);
            % --- Phase loop ----------------------------------------------------
            for ph = 1:nph
                rho_f = op.faceAvg(ifcell(rho, ph));
                s = model.getProp(state, ['s', nm(ph)]);
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
                    v_pore = unorm ./ phi_f;            % interstitial velocity
                    u_dot_n = 0;
                    for d = 1:dim
                        u_dot_n = u_dot_n + (u_face{d} ./ unorm) .* n_unit(:, d);
                    end
                    % n^T D_disp n
                    D_disp = (aT(ph) .* v_pore) + (aL(ph) - aT(ph)) .* v_pore .* (u_dot_n.^2);
                end
                % >>> DIFFUSION PART <<<
                if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                    phiS = phi .* s;
                    tau_MQ = s.*(phiS).^(sf.tortuosityExponent) .* (max(phi, sf.minPorosity)).^(-2);
                else
                    tau_MQ = [];
                end
                % --- Component loop --------------------------------------------
                for c = 1:ncomp
                    z = pick_mole_frac(ph, L_ix, V_ix, xc{c}, yc{c});
                    grad_z_n = op.Grad(z) ./ dist_ij;   % Δz/Δx
                    D_n = D_disp;
                    if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                        if ph == L_ix
                            D_i = tau_MQ .* Dliq_ref(c);          % no rho
                        else
                            invDi = 0;
                            yAll = model.getProps(state, 'y');
                            for j = 1:ncomp
                                if j == c, continue; end
                                yj   = localGetComponentVector(yAll, j);
                                Dij  = Dij_ref(c, j) .* gasScale;   % T & P scaling applied here
                                Dij  = max(Dij, sf.minDiffusivity);
                                invDi = invDi + yj ./ Dij;
                            end
                            % --- Fixed: Standard Wilke Multi-component Numerator (1 - yc) ---
                            yc_curr = localGetComponentVector(yAll, c);
                            Di_mix = max(1 - yc_curr, 1e-6) ./ max(invDi, sf.minDiffusivity);
                            D_i = tau_MQ .* Di_mix;                % no rho
                        end
                        D_i = max(D_i, sf.minDiffusivity);
                        D_n = D_n + op.faceAvg(D_i);
                    end
                    % Final face flux – rho applied once
                    J{c, ph} = - rho_f .* phi_f .* A_int .* D_n .* grad_z_n;
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
function v = localGetComponentVector(z, j)
    if iscell(z), v = z{j}; else v = z(:, j); end
end
function [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(sf, model)
    namecp = model.compFluid.names();
    ncomp  = numel(namecp);
    Dliq_ref = sf.defaultLiquidDiffusivity * ones(ncomp, 1);
    paramLJ  = repmat([sf.defaultSigma, sf.defaultEpsilon], ncomp, 1);
    db.Dliq = struct('h2',6.44e-9,'c1',2.15e-9,'methane',2.15e-9,'co2',2.72e-9, ...
        'h2o',3.29e-9,'water',3.29e-9,'n2',2.86e-9,'c2',1.72e-9, ...
        'ethane',1.72e-9,'c3',1.43e-9,'propane',1.43e-9,'h2s',2.15e-9, ...
        'nc4',1.15e-9,'butane',1.15e-9);
    db.LJ = struct('h2',[2.92,59.7],'c1',[3.758,148.6],'methane',[3.758,148.6], ...
        'co2',[3.996,195.2],'h2o',[2.641,809.1],'water',[2.641,809.1], ...
        'n2',[3.798,71.4],'c2',[4.443,215.7],'ethane',[4.443,215.7], ...
        'c3',[5.118,237.1],'propane',[5.118,237.1],'h2s',[3.60,301.0], ...
        'nc4',[5.206,289.5],'butane',[5.206,289.5]);
    for i = 1:ncomp
        key = lower(namecp{i});
        if isfield(db.Dliq, key), Dliq_ref(i) = db.Dliq.(key); end
        if isfield(db.LJ, key),   paramLJ(i,:) = db.LJ.(key);  end
    end
end
function Dij_ref = localBinaryDiffusionReference(sf, model, paramLJ)
    % Computes the geometric part of binary gas diffusion coefficients
    % scaled directly into SI units (m²/s) at 1 atm.
    ncomp = size(paramLJ, 1);
    Molmass = 1e3 .* model.compFluid.molarMass;   % kg/mol -> g/mol
    sigma = paramLJ(:, 1);                         % collision diameter [Å]
    
    % --- Fixed: Changed from 1.858e-3 (cm²/s) to 1.858e-7 (m²/s) ---
    const = 1.858e-7;
    Dij_ref = zeros(ncomp);
    fTref = sf.Tref^1.5;   % temperature factor at reference conditions
    for i = 1:ncomp
        for j = 1:ncomp
            if i == j, Dij_ref(i,j) = inf; continue; end
            sigma_ij = 0.5 * (sigma(i) + sigma(j));
            sigma_ij2 = sigma_ij^2;
            
            % --- Fixed: Correct Chapman-Enskog mass relationship: sqrt(1/Mi + 1/Mj) ---
            sqrtMass = sqrt((1.0 / Molmass(i)) + (1.0 / Molmass(j)));
            Dij_ref(i,j) = const * fTref*sqrtMass / sigma_ij2;
        end
    end
    Dij_ref = 0.5 * (Dij_ref + Dij_ref.');
end