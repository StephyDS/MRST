%% H2-H2O Vapor-Liquid Equilibrium Calculations
% This script calculates the vapor-liquid equilibrium (VLE) for a hydrogen-water
% mixture under varying conditions using a thermodynamic equation of state.

clear; clc;

% Add necessary MRST modules
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Define Composition Mixture
% Create a compositional mixture with water and hydrogen components
compFluid = TableCompositionalMixture({'Water', 'Hydrogen'}, {'H2O', 'H2'});
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
    case 1 %sources: Solubility of H2 in water and NaCl brine under subsurface 
        % storage conditions (https://hal.science/hal-04623907v1, 2023)
        eosModelsw.msalt=0;
        pres=[100.01, 150.01, 200.01, 200.01, 101.11, 101.31, 130.01,...
            165.01, 199.91, 200.11, 100.01, 100.01, 100.01, 175.11]*barsa;     
        Temp=[298.20,298.05,298.15,298.15,323.55,323.50,323.85,323.55,...
            323.30,323.35,373.85,373.80,373.85,373.65]*Kelvin;
        xliqH2Exp=[0.00135994,0.00199679,0.00264397,0.00263320,0.00125091,...
            0.00123925,0.00159253,0.00201117,0.00244186,0.00245479,...
            0.00142471,0.00140332,0.00140893,0.00243347];

    case 2
        eosModelsw.msalt=1;
        Temp=[298.20,298.30,298.15,298.30,323.20,323.40,323.35,323.20,...
            323.30,323.30,323.20,373.25,373.40,373.10,373.15,373.45,373.15,373.00];
        pres=[100.71,150.01,150.01,200.01,100.31,100.61,101.01,150.06,...
            175.01,199.91,200.01,100.11,126.01,150.36,150.46,150.71,175.51,200.46]*barsa;
        xliqH2Exp=[0.00107012,0.00159298,0.00161244,0.00213575,0.00102099,...
            0.00102078,0.00102426,0.00151377,0.00176728,0.00202590,...
            0.00204487,0.00119671,0.00148492,0.00175595,0.00179573,...
            0.00177350,0.00204000,0.00234604];
    case 3
        eosModelsw.msalt=2;
        Temp=[298.15,298.05,298.05,323.20,323.40,323.35,323.40,373.05,373.20,373.40];
        pres=[100.01,150.01,200.01,100.01,150.01,150.01,200.01,100.01,150.01,200.51]*barsa;
        xliqH2Exp=[0.00088640,0.00132235,0.00171848,0.00088260,0.00131242,...
            0.00128912,0.00172402, 0.00099379, 0.00151866, 0.00205031 ];
    case 4
        eosModelsw.msalt=4;
        Temp=[298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15];
        pres=[100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01]*barsa;
        xliqH2Exp=[0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
            0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];
end

%% case 2, msalt=1
eosModelsw2.msalt=1;
Temp2=[298.20,298.30,298.15,298.30,323.20,323.40,323.35,323.20,...
    323.30,323.30,323.20,373.25,373.40,373.10,373.15,373.45,373.15,373.00];
pres2=[100.71,150.01,150.01,200.01,100.31,100.61,101.01,150.06,...
    175.01,199.91,200.01,100.11,126.01,150.36,150.46,150.71,175.51,200.46]*barsa;
xliqH2Exp2=[0.00107012,0.00159298,0.00161244,0.00213575,0.00102099,...
    0.00102078,0.00102426,0.00151377,0.00176728,0.00202590,...
    0.00204487,0.00119671,0.00148492,0.00175595,0.00179573,...
    0.00177350,0.00204000,0.00234604];

%% case 4, msalt=4
 eosModelsw4 .msalt=4;
 Temp4=[298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15];
 pres4=[100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01]*barsa;
 xliqH2Exp4=[0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
     0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];




%% Perform Flash Calculations
% Determine the liquid-phase hydrogen fraction (xliqH2) for each condition
nc = numel(pres);
nc2 = numel(pres2);
nc4 = numel(pres4);
namecp = eosModelsw.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indH2O= find(strcmp(namecp,'H2O'));

[Lsw, xsw, ~] = standaloneFlash(pres, Temp, z0, eosModelsw); 
[Lsw2, xsw2, ~] = standaloneFlash(pres2, Temp2, z0, eosModelsw2); 
[Lsw4, xsw4, ~] = standaloneFlash(pres4, Temp4, z0, eosModelsw4); 
[Lpr, xpr, ~] = standaloneFlash(pres, Temp, z0, eosModelpr); 
xliqH2sw=xsw(:,indH2);
xliqH2sw2=xsw2(:,indH2);
xliqH2sw4=xsw4(:,indH2);
xliqH2pr=xpr(:,indH2);
presbar=pres./patm;
presbar2=pres2./patm;
presbar4=pres4./patm;

min_xliqH2=[min(xliqH2Exp),min(xliqH2sw),min(xliqH2pr)];
max_xliqH2=[max(xliqH2Exp),max(xliqH2sw),max(xliqH2pr)];

min_xliqH2b=[min(xliqH2Exp2),min(xliqH2Exp4),min(xliqH2sw2),min(xliqH2sw4)];
max_xliqH2b=[max(xliqH2Exp2),max(xliqH2Exp4),max(xliqH2sw2),max(xliqH2sw4)];
minpresbarb=[min(presbar2),min(presbar4)];
maxpresbarb=[max(presbar2),max(presbar4)];

%% calculate the errors
error_swexp=abs(xliqH2sw-xliqH2Exp')./xliqH2Exp';
error_prexp=abs(xliqH2pr-xliqH2Exp')./xliqH2Exp';
fprintf('Errormax SW-Experiment: %12.8f, Errormax PR-Experiment: %12.8f\n', max(error_swexp), max(error_prexp));
fprintf('Errormin SW-Experiment: %12.8f, Errormin PR-Experiment: %12.8f\n', min(error_swexp), min(error_prexp));
fprintf('Errormean SW-Experiment: %12.8f, Errormean PR-Experiment: %12.8f\n', mean(error_swexp), mean(error_prexp));

%% plot the results
f1=figure('Name','H2_solubility_Swmsalt','NumberTitle','off');
f1.Position(3:4) = [900 700];

plot(presbar2,xliqH2sw2,'b*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar4,xliqH2sw4,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar2,xliqH2Exp2,'k o','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar4,xliqH2Exp4,'g diamond','MarkerSize',8,'LineWidth',2)

title('H2 solubility in salt water with the SW model','FontSize',14,'FontWeight','bold','Color','k')
%xlabel({'pressure','(bar)'},'FontSize',14,'FontWeight','bold','Color','k')
%ylabel('H2 molar fraction','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'pressure','(bar)'},'FontWeight','bold','Color','k')
ylabel('H2 molar fraction','FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'x_{H_{2}}, m_{salt}=1','x_{H_{2}, m_{salt}=4}','x_{H_{2}}^{exp}, m_{salt}=1','x_{H_{2}}^{exp}, m_{salt}=4'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')
xlim([min(minpresbarb)-10 max(maxpresbarb)+10])
ylim([min(min_xliqH2b)-1e-4 max(max_xliqH2b)+1e-4])


f2=figure('Name','H2_solubility_SwPr','NumberTitle','off');
f2.Position(3:4) = [900 700];

plot(presbar,xliqH2sw,'b*','MarkerSize',7,'LineWidth',2)
hold on
plot(presbar,xliqH2pr,'r square','MarkerSize',8,'LineWidth',2)
hold on
plot(presbar,xliqH2Exp,'k o','MarkerSize',8,'LineWidth',2)

title('H2 solubility in pure water','FontSize',14,'FontWeight','bold','Color','k')
xlabel({'pressure','(bar)'},'FontWeight','bold','Color','k')
ylabel('H2 molar fraction','FontWeight','bold','Color','k')
ax = gca;
ax.FontSize = 14; 

legend({'x_{H_{2}}, SW','x_{H_{2}}, PR','x_{H_{2}}^{exp}'},...
    'FontSize',13,'TextColor','black',...
    'Location','best')
xlim([min(presbar)-10 max(presbar)+10])
ylim([min(min_xliqH2)-1e-4 max(max_xliqH2)+1e-4])







%% Write Results to File
% Save the results (temperature, pressure, hydrogen mole fraction) to a file
outputFileName = sprintf('H2solubility_case%d_msalt%d_SwPr.dat', caseTest, eosModelsw.msalt);
fileID = fopen(outputFileName, 'wt');

fprintf(fileID, ['Temperature (K)  Pressure (bar)' ...
    ' H2_molar_fraction(SW) H2_molar_fraction(PR)  H2_molar_fraction_exp\n']);
for i = 1:nc
    fprintf(fileID, '%12.2f  %12.4f  %12.8f %12.8f  %12.8f\n',...
        Temp(i), pres(i) / patm, xliqH2sw(i), xliqH2pr(i), xliqH2Exp(i));
end
fclose(fileID);

disp(['Results written to file: ', outputFileName]);


