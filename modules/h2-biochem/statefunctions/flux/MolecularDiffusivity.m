classdef MolecularDiffusivity < StateFunction
    % Compute the effective molecular diffusivity.
    %
    % SYNOPSIS:
    %   d = MolecularDiffusivity(model)
    %   d = MolecularDiffusivity(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Evaluates the effective cell‑centred molecular diffusivity [m²/s]
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance, must have fields `bactDiffusion` (logical)
    %           and `` (positive numeric).
    %
    % OPTIONAL PARAMETERS (property/value pairs):
    %   (none defined in this base class; can be extended via `merge_options`)
    %
    % RETURNS:
    %   dN - Cell‑centred effective diffusivity (ncells × 1), units m²/s.
    %
    % SEE ALSO:
    %   StateFunction, DynamicFlowTransmissibility

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
        %-----------------------------------------------------------------%
        function sf = MolecularDiffusivity(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf = merge_options(sf, varargin{:});

            % Dependencies: pressure, saturation, and nbact (for dynamic porosity)
            sf = sf.dependsOn({'s','pressure'}, 'state');
            sf = sf.dependsOn('temperature', 'state');
            % Dynamic porosity (e.g., bio‑clogging)
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                sf = sf.dependsOn('nbact', 'state');
            end
            sf.label = 'D_{moldiff}^{eff}';
        end

        %-----------------------------------------------------------------%
        function dN = evaluateOnDomain(prop, model, state)
            if ~isfield(state, 'x')
                [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());
                dN = repmat({0}, ncomp, nph);
                return;
            end
            
            [ncomp, nph] = deal(model.getNumberOfComponents(), model.getNumberOfPhases());
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            nm = model.getPhaseNames();   % for saturation look‑up

            % Fetch primary state variables
            [T,p] = model.getProps(state, 'temperature', 'pressure');

            % Porosity – dynamic (bio‑clogging) or static
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                nbact = model.getProp(state, 'nbact');
                phi = model.rock.poro(p, nbact);
            else
                phi = model.rock.poro;
            end
            % --- Diffusion parameters (if enabled) ------------------------------
            p_safe = max(p, 1e-8*barsa);              % ensure pressure > 0
            gasScale = (T./prop.Tref).^prop.gasDiffExponent .* (prop.pref./p_safe);
            [Dliq_ref, paramLJ] = localLoadDiffusionDatabase(prop, model);
            Dij_ref = localBinaryDiffusionReference(prop, model, paramLJ);  % no temperature factor inside

           
            

            % --- Output initialisation ------------------------------------------
            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                dN = cell(ncomp, nph); [dN{:}] = deal(0);
                for ph = 1:nph
                    s = model.getProp(state, ['s', nm(ph)]);
                    Voln = max(s, 1.0e-8);
                    phiS = phi .* s;
                    tau_MQ = (phiS).^(prop.tortuosityExponent) .* (max(phi, prop.minPorosity)).^(-2);

                    % --- Component loop --------------------------------------------
                    for c = 1:ncomp
                        if ph == L_ix
                            D_i = tau_MQ .* Dliq_ref(c);          % no rho
                        elseif ph == V_ix
                            invDi = 0;
                            yAll = model.getProps(state, 'y');
                            for j = 1:ncomp
                                if j == c, continue; end
                                yj   = localGetComponentVector(yAll, j);
                                Dij  = Dij_ref(c, j) .* gasScale;   % T & P scaling applied here
                                Dij  = max(Dij, prop.minDiffusivity);
                                invDi = invDi + yj ./ Dij;
                            end
                            % --- Fixed: Standard Wilke Multi-component Numerator (1 - yc) ---
                            yc_curr = localGetComponentVector(yAll, c);
                            Di_mix = max(1 - yc_curr, 1e-6) ./ max(invDi, prop.minDiffusivity);
                            D_i = tau_MQ .* Di_mix;                % no rho
                        end
                        D_i = max(D_i, prop.minDiffusivity);

                        % Final face flux – rho applied once
                        dN{c, ph} = phi .*Voln .* D_i;
                    end
                end
            else
                [dN{:}] = deal(0);
            end
        end
    end
end


% ------------------------------------------------------------------------
% Helper functions
% ------------------------------------------------------------------------
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