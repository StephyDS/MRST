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
            dispEffects = (isprop(model, 'molecularDispersion') && model.molecularDispersion)||...
                (isprop(model, 'molecularDiffusion') && model.molecularDiffusion);
            if dispEffects
                sf = sf.dependsOn('ComponentTotalMolecularDispFlux');
            end
            sf.label = 'V_i';
        end

        function v = evaluateOnDomain(sf, model, state)
            % 1. Get the standard advective total flux (cell array per component)
            v = evaluateOnDomain@ComponentTotalFlux(sf, model, state);

            ncomp = numel(v);   % number of components

            % 2. Optionally add mechanical dispersion or diffusion
            dispEffects = (isprop(model, 'molecularDispersion') && model.molecularDispersion)||...
                (isprop(model, 'molecularDiffusion') && model.molecularDiffusion);
            if dispEffects
                Jdispdiff = sf.getEvaluatedDependencies(state, 'ComponentTotalMolecularDispFlux');
                for c = 1:ncomp
                    if numel(Jdispdiff) >= c && ~isempty(Jdispdiff{c})
                        v{c} = v{c} + Jdispdiff{c};
                    end
                end
            end
        end
    end
end