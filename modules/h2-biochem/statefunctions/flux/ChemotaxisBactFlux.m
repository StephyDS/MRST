classdef ChemotaxisBactFlux < StateFunction
    % Chemotaxis flux of bacteria (mass‑fraction based).
    %
    % SYNOPSIS:
    %   flux = ChemotaxisBactFlux(model)
    %
    % DESCRIPTION:
    %   Computes the face‑wise Chemotaxis flux of bacteria in the liquid
    %   phase as
    % Jchemotaxis = Kchemo * nbact *(1 - nbact) * Grad(xH2)
    % (volume-filling effect, generalized Keller-Segel model)
   %   where T is the chemotaxis transmissibility provided by
    %   `chemotaxisTransmissibility` and ρ_l^f is the face‑averaged liquid
    %   density.
    %
    % REQUIRED PARAMETERS:
    %   model - Model with `chemotaxisTransmissibility` state function
    %           registered and bacterial modelling enabled.
    %
    % SEE ALSO:
    %   StateFunction, chemotaxisTransmissibility

    properties
    end

    methods
        %-----------------------------------------------------------------%
        function df = ChemotaxisBactFlux(model, varargin)
            df@StateFunction(model, varargin{:});
            df = merge_options(df, varargin{:});

            % Dependencies
            df = df.dependsOn({'ChemotaxisTransmissibility'});
            df = df.dependsOn({'x'}, 'state');
            df = df.dependsOn('Density', 'PVTPropertyFunctions');

            df.label = 'J_{bact}^{chemo}';
        end

        %-----------------------------------------------------------------%
        function DN = evaluateOnDomain(prop, model, state)
            op = model.operators;
            bcrm=model.biochemFluid;

            % Retrieve microbial transmissibility
            T = prop.getEvaluatedDependencies(state, 'ChemotaxisTransmissibility');

            % State variables
            [x, rho] = model.getProps(state, 'x', 'Density');

            % Liquid density (face averaged)
            namecp = model.getComponentNames();
            indH2= find(strcmpi(namecp, bcrm.rH2), 1);
            L_ix = model.getLiquidIndex();
            if iscell(rho)
                rhoL = rho{L_ix};
            else
                rhoL = rho(:, L_ix);
            end
            rhoLf = model.operators.faceAvg(rhoL);

            if iscell(x)
                xH2=x{indH2};
            else
                xH2=x(:,indH2);
            end

            % Diffusive flux (all faces)
            DN = + rhoLf .* T .* op.Grad(xH2);
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