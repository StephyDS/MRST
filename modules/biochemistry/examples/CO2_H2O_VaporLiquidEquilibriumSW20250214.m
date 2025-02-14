%% CO2-H2O Vapor-Liquid Equilibrium Calculations
% This script calculates the vapor-liquid equilibrium (VLE) for a hydrogen-water
% mixture under varying conditions using a thermodynamic equation of state.

clear; clc;

% Add necessary MRST modules
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'CarbonDioxide'}, {'H2O', 'CO2'});
disp(compFluid);

%% Initialize Thermodynamic Model
% Choose an equation of state model for the calculations
eosNamesw = 'sw'; % Soreide-Whitson (SW) model
eosModelsw = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamesw);
eosModelsw2 = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamesw);
eosModelsw4 = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamesw);

eosNamepr = 'pr'; % Soreide-Whitson (SW) model
eosModelpr = SoreideWhitsonEquationOfStateModel([], compFluid, eosNamepr);

%% Define Test Case Parameters
% Set the test case for different pressures, temperatures, and salinity levels
z0 = [0.8, 0.2]; % Initial composition
patm = 1e5; % Atmospheric pressure in Pa
caseTest = 1; % Choose the test case here

switch caseTest
    case 1 %sources: Solubility of CO2 in water and NaCl brine under subsurface 
        % storage conditions (https://hal.science/hal-04623907v1, 2023)
        eosModelsw.msalt=0;
        Temp=[323.2, 323.2, 323.2, 323.2, 323.2, 323.2, 333.2, 333.2, 333.2,...
            333.2, 333.2, 333.2, 353.1, 353.1, 353.1, 353.1, 353.1, 353.1]*Kelvin;
        pres=[40.5, 60.6, 80.8, 100.9, 121, 141.1, 40.5, 60.6, 80.8, 100.9,...
            121, 141.1, 40.5, 60.6, 80.8, 100.9, 121, 131]*barsa;
        xliqCO2Exp=[0.0109, 0.0161, 0.019, 0.0205, 0.0214, 0.0217, 0.0096,...
            0.0138, 0.0166, 0.0186, 0.0201, 0.0208, 0.008, 0.0114, 0.014, ...
            0.016, 0.0176, 0.0184];

    case 2
        eosModelsw.msalt=1;
        Temp=[373.38, 373.37, 373.41]*Kelvin;
        pres=[16.983, 32.527, 68.182]*barsa;
        xliqCO2Exp=[0.00237, 0.00426, 0.00833];
    case 3
         eosModelsw.msalt=1.13;
        Temp=[323.02, 322.97, 323.03, 323.04, 372.33, 372.31, 372.29, 372.29, 372.25]*Kelvin;
        pres=[53.450, 75.550, 100.350, 145.080, 31.148, 60.500, 108.840, 151.920, 191.980]*barsa;
        xliqCO2Exp=[0.01030, 0.01290, 0.01510, 0.01700,0.00390, 0.00750, 0.01130, 0.01360, 0.01570];
            
    case 4
        eosModelsw.msalt=3.01;
        Temp=[342.82, 342.81, 342.82, 372.39, 372.42, 372.41 , 372.43, 372.45, 372.45]*Kelvin;
        pres=[30.391, 72.559, 100.910, 25.556, 71.417, 100.517, 152.433, 199.597, 229.817]*barsa;
        xliqCO2Exp=[0.00441, 0.00880, 0.01057, 0.00292, 0.00707, 0.00878, 0.01141, 0.01258, 0.01337];
end

%% case 2, msalt=1
eosModelsw2.msalt=1;
Temp2=[373.38, 373.37, 373.41]*Kelvin;
pres2=[16.983, 32.527, 68.182]*barsa;
xliqCO2Exp2=[0.00237, 0.00426, 0.00833];

%% case 4, msalt=4
  eosModelsw4.msalt=3.01;
  Temp4=[342.82, 342.81, 342.82, 372.39, 372.42, 372.41 , 372.43, 372.45, 372.45]*Kelvin;
  pres4=[30.391, 72.559, 100.910, 25.556, 71.417, 100.517, 152.433, 199.597, 229.817]*barsa;
  xliqCO2Exp4=[0.00441, 0.00880, 0.01057, 0.00292, 0.00707, 0.00878, 0.01141, 0.01258, 0.01337];



%% Perform Flash Calculations
% Determine the liquid-phase hydrogen fraction (xliqCO2) for each condition
nc = numel(pres);
nc2 = numel(pres2);
nc4 = numel(pres4);
namecp = eosModelsw.getComponentNames();
indCO2=find(strcmp(namecp,'CO2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pres, Temp, z0, eosModelsw); 
[Lsw2, xsw2, ~] = standaloneFlash(pres2, Temp2, z0, eosModelsw2); 
[Lsw4, xsw4, ~] = standaloneFlash(pres4, Temp4, z0, eosModelsw4); 
[Lpr, xpr, ~] = standaloneFlash(pres, Temp, z0, eosModelpr); 
xliqCO2sw=xsw(:,indCO2);
xliqCO2sw2=xsw2(:,indCO2);
xliqCO2sw4=xsw4(:,indCO2);
xliqCO2pr=xpr(:,indCO2);
presbar=pres./patm;
presbar2=pres2./patm;
presbar4=pres4./patm;

min_xliqCO2=[min(xliqCO2Exp),min(xliqCO2sw),min(xliqCO2pr)];
max_xliqCO2=[max(xliqCO2Exp),max(xliqCO2sw),max(xliqCO2pr)];

min_xliqCO2b=[min(xliqCO2Exp2),min(xliqCO2Exp4),min(xliqCO2sw2),min(xliqCO2sw4)];
max_xliqCO2b=[max(xliqCO2Exp2),max(xliqCO2Exp4),max(xliqCO2sw2),max(xliqCO2sw4)];
minpresbarb=[min(presbar2),min(presbar4)];
maxpresbarb=[max(presbar2),max(presbar4)];

%% calculate the errors
error_swexp=abs(xliqCO2sw-xliqCO2Exp')./xliqCO2Exp';
error_prexp=abs(xliqCO2pr-xliqCO2Exp')./xliqCO2Exp';
fprintf('Errormax SW-Experiment: %12.8f, Errormax PR-Experiment: %12.8f\n', max(error_swexp), max(error_prexp));
fprintf('Errormin SW-Experiment: %12.8f, Errormin PR-Experiment: %12.8f\n', min(error_swexp), min(error_prexp));
fprintf('Errormean SW-Experiment: %12.8f, Errormean PR-Experiment: %12.8f\n', mean(error_swexp), mean(error_prexp));

%% plot the results
f1=figure('Name','CO2_solubility_Swmsalt','NumberTitle','off');
f1.Position(3:4) = [900 700];

plot(presbar2,xliqCO2sw2,'b*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar4,xliqCO2sw4,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar2,xliqCO2Exp2,'k o','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar4,xliqCO2Exp4,'g diamond','MarkerSize',8,'LineWidth',2)

title('CO2 solubility in salt water with the SW model','FontSize',14,'FontWeight','bold','Color','k')
%xlabel({'pressure','(bar)'},'FontSize',14,'FontWeight','bold','Color','k')
%ylabel('CO2 molar fraction','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'pressure','(bar)'},'FontWeight','bold','Color','k')
ylabel('CO2 molar fraction','FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'x_{CO_{2}}, m_{salt}=1','x_{CO_{2}, m_{salt}=3.01}','x_{CO_{2}}^{exp}, m_{salt}=1','x_{CO_{2}}^{exp}, m_{salt}=3.01'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')
xlim([min(minpresbarb)-10 max(maxpresbarb)+10])
ylim([min(min_xliqCO2b)-1e-4 max(max_xliqCO2b)+1e-4])


f2=figure('Name','CO2_solubility_SwPr','NumberTitle','off');
f2.Position(3:4) = [900 700];

plot(presbar,xliqCO2sw,'b*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar,xliqCO2pr,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar,xliqCO2Exp,'k o','MarkerSize',8,'LineWidth',2)

title('CO2 solubility in pure water','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'pressure','(bar)'},'FontWeight','bold','Color','k')
ylabel('CO2 molar fraction','FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'x_{CO_{2}}, SW','x_{CO_{2}}, PR','x_{CO_{2}}^{exp}'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')
xlim([min(presbar)-10 max(presbar)+10])
ylim([min(min_xliqCO2)-1e-4 max(max_xliqCO2)+1e-4])







%% Write Results to File
% Save the results (temperature, pressure, hydrogen mole fraction) to a file
outputFileName = sprintf('CO2solubility_case%d_msalt%d_SwPr.dat', caseTest, eosModelsw.msalt);
fileID = fopen(outputFileName, 'wt');

fprintf(fileID, ['Temperature (K)  Pressure (bar)' ...
    ' CO2_molar_fraction(SW) CO2_molar_fraction(PR)  CO2_molar_fraction_exp\n']);
for i = 1:nc
    fprintf(fileID, '%12.2f  %12.4f  %12.8f %12.8f  %12.8f\n',...
        Temp(i), pres(i) / patm, xliqCO2sw(i), xliqCO2pr(i), xliqCO2Exp(i));
end
fclose(fileID);

disp(['Results written to file: ', outputFileName]);


