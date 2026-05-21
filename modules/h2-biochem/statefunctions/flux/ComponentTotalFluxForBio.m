classdef ComponentTotalFluxForBio < ComponentTotalFlux
    % Total component flux for bio‑chemical models.
    % Inherits the  ComponentTotalFlux
    % and optionally adds molecular dispersion and diffusion fluxes
    % when the corresponding model flags are enabled:
    %   model.molecularDispersion  (logical)
    %   model.molecularDiffusion   (logical)

    methods
        function sf = ComponentTotalFluxForBio(model, varargin)
            % Let the parent class register its dependencies (ComponentPhaseFlux, etc.)
            sf = sf@ComponentTotalFlux(model, varargin{:});

            % Add optional dependencies for the extra physics
            if isprop(model, 'molecularDispersion') && model.molecularDispersion
                sf = sf.dependsOn('ComponentTotalMolecularDispFlux');
            end
            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                sf = sf.dependsOn('ComponentTotalMolecularDiffFlux');
            end

            sf.label = 'V_i';
        end

        function v = evaluateOnDomain(sf, model, state)
            % 1. Get the standard advective total flux (cell array per component)
            v = evaluateOnDomain@ComponentTotalFlux(sf, model, state);

            ncomp = numel(v);   % number of components

            % 2. Optionally add mechanical dispersion
            if isprop(model, 'molecularDispersion') && model.molecularDispersion
                Jdisp = sf.getEvaluatedDependencies(state, 'ComponentTotalMolecularDispFlux');
                for c = 1:ncomp
                    if numel(Jdisp) >= c && ~isempty(Jdisp{c})
                        v{c} = v{c} + Jdisp{c};
                    end
                end
            end

            % 3. Optionally add molecular diffusion
            if isprop(model, 'molecularDiffusion') && model.molecularDiffusion
                Jdiff = sf.getEvaluatedDependencies(state, 'ComponentTotalMolecularDiffFlux');
                for c = 1:ncomp
                    if numel(Jdiff) >= c && ~isempty(Jdiff{c})
                        v{c} = v{c} + Jdiff{c};
                    end
                end
            end
        end
    end
end