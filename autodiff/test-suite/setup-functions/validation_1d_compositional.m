function setup = validation_1d_compositional(varargin)
% Setup function for a quarter five-spot water-oil model
%
% SYNOPSIS:
%   setup = validation_1d_compositional('pn1', pv1, ...)
%   setup = validation_1d_compositional(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a validation case from the compositional module that studies
%   the injection of CO2 into a 1D reservoir model filled with decane, CO2,
%   and methane.
%
%   Configurable parameters: none.
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

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
    % One-line description
    description ...
        = ['Example from the compositional module with CO2 injection ', ...
           'into 1D reservoir filled with Decane, CO2 and Methane'    ];
       
    % Optional input arguments
    options = struct();
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props compositional deckformat
    
    % Read deck    
    pth = getDatasetPath('simplecomp');
    fn  = fullfile(pth, 'SIMPLE_COMP.DATA');
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);

    % Set up grid
    G = initEclipseGrid(deck);
    G = computeGeometry(G);
    
    % Set up rock
    rock  = initEclipseRock(deck);
    rock  = compressRock(rock, G.cells.indexMap);
    fluid = initDeckADIFluid(deck);
    
    % Define some surface densities
    fluid.rhoOS = 800;
    fluid.rhoGS = 10;
    
    % Define model
    eos = initDeckEOSModel(deck);
    model = GenericOverallCompositionModel(G, rock, fluid, eos, 'water', false);
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', false);
    
    % Get schedule
    schedule = convertDeckScheduleToMRST(model, deck);
    
    % Manually set the injection composition
    [schedule.control.W.components] = deal([0, 1, 0]);
    
    % Injection is pure gas
    [schedule.control.W.compi] = deal([1, 0]);
    
    % The problem is defined at 150 degrees celsius with 75 bar initial
    % pressure. We set up the initial problem and make a call to the flash
    % routines to get correct initial composition.
    for i = 1:numel(schedule.control.W)
        schedule.control.W(i).lims = [];
    end
    
    % Initial conditions
    z0 = [0.6, 0.1, 0.3];
    T  = 150 + 273.15;
    p  = 75*barsa;
    state0 = initCompositionalState(G, p, T, [1, 0], z0, eos);
    
    % Plotting
    plotOptions = {'plot1d'            , true          , ...
                   'lockCaxis'         , true          , ...
                   'Size'              , [800, 350]    , ...
                   'PlotBoxAspectRatio', [2.5,1,1]     , ...
                   'Projection'        , 'Orthographic', ...
                   'YLim'              , [0,1]         , ...
                   'View'              , [0,90]        };
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end