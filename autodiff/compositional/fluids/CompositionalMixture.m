classdef CompositionalMixture
    properties
        Tcrit % Critical temperature in Kelvin
        Pcrit % Critical pressure in Pascal
        Vcrit % Critical volume in m^3 / mol
        acentricFactors % Acentric factors (dimensionless)
        molarMass % Component mass (kg / mol)
        names % Names of each component. Each name must be unique.
        name % Name of the fluid mixture itself
    end
    
    properties ( Access = protected)
        bic % Binary interaction coefficients
    end
    
    methods
        function fluid = CompositionalMixture(names, Tcrit, Pcrit, Vcrit, acentricFactors, molarMass, varargin)
            if nargin == 0
                disp('Empty fluid - no validation');
                return
            end
            
            for i = 1:numel(names)
                cts = sum(strcmp(names{i}, names));
                if cts > 1
                    warning(['Component ', num2str(i), ': ', names{i}, ' occurs multiple times.']);
                end
            end
            fluid.names = names;
            fluid.Tcrit = Tcrit;
            fluid.Pcrit = Pcrit;
            fluid.Vcrit = Vcrit;
            fluid.acentricFactors = acentricFactors;
            fluid.molarMass = molarMass;
            ncomp = numel(fluid.names);
            
            fluid.bic = zeros(ncomp, ncomp);
        end
        
        function n = getNumberOfComponents(fluid)
            n = numel(fluid.names);
        end
        
        function bic = getBinaryInteraction(fluid)
            bic = fluid.bic;
        end
  
%=============SDS MODIF=================================
        function bic = getBinaryInteractionGH2H2O(fluid,T)
                       
            namecp=fluid.names;
            indH2O=find(strcmp(namecp,'H2O'));
            indH2=find(strcmp(namecp,'H2'));
            indCO2=find(strcmp(namecp,'CO2'));
            
            if ~isempty(indH2)
                TrH2=T./fluid.Tcrit(indH2);
                [D0,D1] = deal(0.01993,0.042834);
                bic{indH2,indH2O} = D0+D1.*TrH2;
                bic{indH2O,indH2} = D0+D1.*TrH2;
                bic{indH2,indH2} = 0.*TrH2 + fluid.bic(indH2,indH2);
                bic{indH2O,indH2O} =0.*TrH2 +fluid.bic(indH2O,indH2O);
            end
            if ~isempty(indCO2)
                TrCO2=T./fluid.Tcrit(indCO2);
                [D0,D1] = deal(0.00068208385571,0.02066623464504);
                bic{indCO2,indH2O} = D0.*TrCO2-D1;
                bic{indH2O,indCO2} = D0.*TrCO2-D1;
                bic{indCO2,indCO2} = 0.*TrCO2 + fluid.bic(indCO2,indCO2);
                bic{indH2O,indH2O} =0.*TrCO2 +fluid.bic(indH2O,indH2O);
            end

        end

        % function bic = getBinaryInteractionGH2H2O(fluid,T)
        % 
        %     namecp=fluid.names;
        %     indH2=find(strcmp(namecp,'H2'));
        %     indH2O=find(strcmp(namecp,'H2O'));
        %     TrH2=T./fluid.Tcrit(indH2);
        %     [D0,D1] = deal(0.01993,0.042834);
        %     bic{indH2,indH2O} = D0+D1.*TrH2;
        %     bic{indH2O,indH2} = D0+D1.*TrH2;
        %     bic{indH2,indH2} = 0.*TrH2 + fluid.bic(indH2,indH2);
        %     bic{indH2O,indH2O} = 0.*TrH2 + fluid.bic(indH2O,indH2O);
        % end

        function bic = getBinaryInteractionLH2H2O(fluid,T,msalt)
            namecp=fluid.names;            
            indH2O=find(strcmp(namecp,'H2O'));
            indH2=find(strcmp(namecp,'H2'));            
            indCO2=find(strcmp(namecp,'CO2'));

            bic{indH2O,indH2O} = fluid.bic(indH2O,indH2O);
            if ~isempty(indH2)
                TrH2=T./fluid.Tcrit(indH2);
                [D0,D1,D2,D3] = deal(-2.11917,0.14888,-13.01835,-0.43946);
                [a0,a1] = deal(-2.26322e-2,-4.4736e-3);
                bic{indH2,indH2O} = D0.*(1+a0.*msalt)+D1.*TrH2*(1+a1*msalt)+D2*exp(D3.*TrH2);
                bic{indH2O,indH2} = D0.*(1+a0.*msalt)+D1.*TrH2*(1+a1*msalt)+D2*exp(D3.*TrH2);
                bic{indH2,indH2} = 0.*TrH2 + fluid.bic(indH2,indH2);
                bic{indH2O,indH2O} =0.*TrH2 +fluid.bic(indH2O,indH2O);
            end
            if ~isempty(indCO2)
                msalt2=msalt*msalt;
                TrCO2=T./fluid.Tcrit(indCO2);
                [D0,D1,D2] = deal(0.43575155,-0.05766906744,0.00826464849);
                [D3,D4,D5] = deal(0.00129539193,-0.0016698848,-0.47866096);
                bic{indCO2,indH2O} = TrCO2.*(D0+D1.*TrCO2+msalt*D2.*TrCO2)...
                    +msalt2.*(D3+D4.*TrCO2)+D5;
                bic{indH2O,indCO2} = TrCO2.*(D0+D1.*TrCO2+msalt*D2.*TrCO2)...
                    +msalt2.*(D3+D4.*TrCO2)+D5;
                bic{indCO2,indCO2} = 0.*TrCO2 + fluid.bic(indCO2,indCO2);
                bic{indH2O,indH2O} =0.*TrCO2 +fluid.bic(indH2O,indH2O);
            end            
        end
        %===========SDS MODIF====================



        function fluid = setBinaryInteraction(fluid, input)
            % Set BIC via a matrix. Must be symmetric and ncomp by ncomp
            ncomp = fluid.getNumberOfComponents();
            
            if isvector(input)
                n_el = ncomp*(ncomp-1)/2;
                assert(numel(input) == n_el);
                
                tmp = zeros(ncomp, ncomp);
                offset = 0;
                for i = 2:ncomp
                    cts = i-1;
                    tmp(i, 1:i-1) = input(offset + (1:cts));
                    offset = offset + cts;
                end
                input = tmp + tmp';
            end
            assert(size(input, 1) == ncomp)
            assert(size(input, 1) == size(input, 2));
            assert(all(all(input == input')));
            fluid.bic = input;
        end
        
        function disp(mixture)
            % builtin('disp', mixture);
            cn = class(mixture);
            isDesktop = usejava('desktop');
            if isDesktop
                fprintf('  <a href="matlab:helpPopup %s">%s</a>:\n', cn, cn);
            else
                fprintf('  %s:\n', cn);
            end
            fprintf('\n');
            nc = mixture.getNumberOfComponents();
            cnames = mixture.names;
            fprintf('  %d component mixture', nc);
            if ~isempty(mixture.name)
                fprintf(' (%s)', mixture.name)
            end
            fprintf(':\n');

            ml = max(max(cellfun(@numel, cnames)), 4);
            sep = repmat('-', 55 + ml, 1);
            fprintf('  %*s | p_c [Pa] | T_c [K] | V_c [m^3] |  acf  | mw [kg/mol] \n', ml, 'Name');
            fprintf('  %s\n', sep)
            for i = 1:nc
                pc = mixture.Pcrit(i);
                tc = mixture.Tcrit(i);
                vc = mixture.Vcrit(i);
                acf = mixture.acentricFactors(i);
                mw = mixture.molarMass(i);
                fprintf('  %*s | %1.2e | %3.1f K | %2.3e | %1.3f | %1.7f \n', ...
                        ml, cnames{i}, pc, tc, vc, acf, mw);
            end
            fprintf('  %s\n', sep);
            bi = mixture.getBinaryInteraction();
            fprintf('  ');
            if any(bi(:))
                fprintf('Binary interaction coefficients:\n');
                disp(bi);
            else
                fprintf('No non-zero binary interaction coefficients.\n')
            end
            fprintf('\n');
            known = {'Pcrit', 'Tcrit', 'Vcrit', 'acentricFactors', 'molarMass', 'names', 'name'};
            props = propertynames(mixture);
            extra = setdiff(props, known);
            if numel(extra)
                fprintf('\n  Additional properties:\n');
                for i = 1:numel(extra)
                    e = extra{i};
                    fprintf('  %s (%s)\n', e, class(mixture.(e)));
                end
                fprintf('\n');
            end
        end
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
