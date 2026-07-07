classdef MicrobialChemotaxis < StateFunction
    % Compute the effective microbial diffusivity in the liquid phase.
    %
    % SYNOPSIS:
    %   d = MicrobialChemotaxis(model)
    %   d = MicrobialChemotaxis(model, 'pn', pv, ...)
    %
    % DESCRIPTION:
    %   Evaluates the effective cell‑centred microbial Chemotaxis [m²/s]
    %
    % REQUIRED PARAMETERS:
    %   model - Model instance, must have fields `bactDiffusion` (logical)
    %           and `biochemFluid.bactdiff` (positive numeric).
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
        % No additional properties needed
    end

    methods
        %-----------------------------------------------------------------%
        function md = MicrobialChemotaxis(model, varargin)
            md@StateFunction(model, varargin{:});
            md = merge_options(md, varargin{:});

            % Dependencies: pressure, saturation, and nbact (for dynamic porosity)
            md = md.dependsOn({'s','nbact'}, 'state');
            %md = md.dependsOn('PoreVolume', 'PVTPropertyFunctions');
            md.label = 'D_{chemot}^{eff}';
        end

        %-----------------------------------------------------------------%
        function dN = evaluateOnDomain(prop, model, state)
            bcrm=model.biochemFluid;
            nbioreact=bcrm.nbioreact;
            % Fetch primary state variables
            [ s, nbact] = model.getProps(state, 's','nbact');

            % Porosity – dynamic (bio‑clogging) or static
            if isprop(model, 'rock') && isa(model.rock.poro, 'function_handle')
                phi = model.rock.poro(p, nbact);
                nbioreact=model.biochemFluid.nbioreact;
                if nbioreact==1
                    if iscell(nbact)
                        phi = model.rock.poro(p, nbact{1}); % Apply both modifications
                    else
                        phi = model.rock.poro(p, nbact);
                    end
                elseif nbioreact==2
                    if iscell(nbact)
                        phi = model.rock.poro(p, nbact{1}, nbact{2}); % Apply both modifications
                    else
                        phi = model.rock.poro(p, nbact(:,1), nbact(:,2)); % Apply both modifications
                    end
                end
            else
                phi = model.rock.poro;
            end


            % Liquid saturation
            L_ix = model.getLiquidIndex();
            if iscell(s)
                sL = s{L_ix};
            else
                sL = s(:, L_ix);
            end

            Voln = max(sL, 1.0e-8);

            % Initialize with zero growth rate
            dN=cell(1,nbioreact);
            [dN{:}] = deal(0);

            % Effective microbial diffusivity
            for i=1:nbioreact
                if iscell(nbact)
                    nbacti=nbact{i};
                else
                    nbacti=nbact(:,i);
                end
                % Effective microbial diffusivity
                if model.chemotaxisEffect && isprop(model.biochemFluid, 'xch_seuil')
                    dN{i} = model.biochemFluid.xch_seuil(i) .* phi .* Voln.*nbacti.*(1.0-nbacti);
                else
                    dN{i} = 0;
                end
            end
        end
    end
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