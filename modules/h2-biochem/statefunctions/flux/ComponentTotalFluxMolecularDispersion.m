classdef ComponentTotalFluxMolecularDispersion < StateFunction
    % Total component flux including advection + molecular diffusion:
    %   v_i = sum_alpha v_{i,alpha} + J_i^{mol,diff}

    methods
        function sf = ComponentTotalFluxMolecularDispersion(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf = sf.dependsOn('ComponentPhaseFlux');
            sf = sf.dependsOn('ComponentTotalMolecularDispFlux');
            sf.label = 'V_i';
        end

        function v = evaluateOnDomain(sf, model, state)
            ncomp = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();

            phaseFlux = sf.getEvaluatedDependencies(state, 'ComponentPhaseFlux');
            Jdisp     = sf.getEvaluatedDependencies(state, 'ComponentTotalMolecularDispFlux');

            v = cell(ncomp, 1);
            for c = 1:ncomp
                vc = 0;
                for ph = 1:nph
                    f = phaseFlux{c, ph};
                    if ~isempty(f)
                        vc = vc + f;
                    end
                end
                if ~isempty(Jdisp) && numel(Jdisp) >= c && ~isempty(Jdisp{c})
                    vc = vc + Jdisp{c};
                end
                v{c} = vc;
            end
        end
    end
end