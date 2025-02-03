function schedule = createCyclicScenario(TotalTime, rampupTime, nCycles, buildUpLength, idleLength, injectionLength, productionLength, W)
    % Create time steps
    deltaT = rampupTimesteps(TotalTime, rampupTime, 0);
    schedule = simpleSchedule(deltaT);
    
    % Initialize control steps
    schedule.step.control = zeros(size(deltaT));
    tmp = cell(4,1);
    schedule.control = struct('W',tmp);
    % Assign control 1 for the build-up phase
    schedule.step.control(1:buildUpLength) = 1;
    schedule.control(1).W = W(1);
    % Assign cyclic controls for n cycles
    currentStep = buildUpLength + 1;
    for cycle = 1:nCycles
        schedule.step.control(currentStep:currentStep+idleLength-1) = 2;
        schedule.control(2).W = W(2);
        currentStep = currentStep + idleLength;
        schedule.step.control(currentStep:currentStep+injectionLength-1) = 3;
        schedule.control(3).W = W(3);
        currentStep = currentStep + injectionLength;
        schedule.step.control(currentStep:currentStep+productionLength-1) = 4;
        schedule.control(4).W = W(4);
        currentStep = currentStep + productionLength;
    end
    
    % Ensure the schedule does not exceed the total number of steps
    schedule.step.control = schedule.step.control(1:length(deltaT));
end