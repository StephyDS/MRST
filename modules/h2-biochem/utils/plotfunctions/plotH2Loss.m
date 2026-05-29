function plotH2Loss(model, schedule, states, ws)
% plotH2Loss  Compute and plot cumulative H2 loss (% of total introduced)
%
%   Uses the well's built-in .H2 field (kg/s) which is:
%       positive  → mass added to reservoir (injection)
%       negative  → mass removed from reservoir (production)
%   Loss = 100 * (initialMass + cumInjected - cumProduced - currentMass)
%               / (initialMass + cumInjected)

% Index of Hydrogen
compNames = model.compFluid.names;
iH2 = find(strcmpi(compNames, 'Hydrogen'));
if isempty(iH2), iH2 = 2; end

nSteps = numel(states);
dt     = schedule.step.val;         % seconds
timeDays = cumsum(dt) / day;

% Current total H2 mass from state
currentMass = zeros(nSteps, 1);
for k = 1:nSteps
    currentMass(k) = sum(states{k}.FlowProps.ComponentTotalMass{iH2});
end

% Cumulative injected and produced H2 mass (using .H2 field)
injectedMassStep = zeros(nSteps, 1);
producedMassStep = zeros(nSteps, 1);

for k = 1:nSteps
    for w = 1:numel(ws{k})
        well = ws{k}(w);
        if isempty(well) || ~isfield(well, 'H2') || isempty(well.H2)
            continue;
        end
        h2Rate = well.H2;   % kg/s (positive = into reservoir)
        if h2Rate > 0
            injectedMassStep(k) = injectedMassStep(k) + h2Rate;
        else
            producedMassStep(k) = producedMassStep(k) - h2Rate;   % make positive
        end
    end
end

cumInjected = cumsum(injectedMassStep .* dt);
cumProduced = cumsum(producedMassStep .* dt);

initialMass    = currentMass(1);
totalIntroduced = initialMass + cumInjected;
consumed        = totalIntroduced - cumProduced - currentMass;   % bacterial consumption
lossPercent     = max(0, 100 * consumed ./ max(totalIntroduced, 1e-12));

% ----- Diagnostic plot -----
figure;
subplot(2,1,1);
plot(timeDays, currentMass, 'b', ...
    timeDays, totalIntroduced - cumProduced, 'r--', 'LineWidth',1.5);
legend('Current mass', 'Expected mass (init+inj-prod)', 'Location','best');
xlabel('Time (days)'); ylabel('H_2 mass (kg)');
title('Mass balance check');
grid on;

subplot(2,1,2);
plot(timeDays, cumInjected, 'g', ...
    timeDays, cumProduced, 'm', ...
    timeDays, cumInjected - cumProduced, 'k', 'LineWidth',1.5);
legend('Cum. injected', 'Cum. produced', 'Net injected', 'Location','best');
xlabel('Time (days)'); ylabel('Cumulative mass (kg)');
grid on;

% ----- Loss percentage plot -----
figure;
plot(timeDays, lossPercent, 'LineWidth',1.5);
xlabel('Time (days)');
ylabel('H_2 loss (%)');
title('Cumulative H_2 loss (percent of total introduced)');
grid on;
end