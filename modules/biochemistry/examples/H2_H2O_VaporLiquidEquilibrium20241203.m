%% H2-H2O Vapor-Liquid Equilibrium Calculations
clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui

%% Composition Mixture
%compFluid =TableCompositionalMixture({'Water','CarbonDioxide'}, {'H2O','CO2'});
compFluid =TableCompositionalMixture({'Water','Hydrogen'}, {'H2O','H2'});
%compFluid = compFluid.setBinaryInteraction(bic);
disp (compFluid)

%% Updating the Thermodynamic Equilibrium
% To use the Rachford-Rice method to determine the vapor-liquid
% equilibrium, we must know the K-values and suitable relationships for the
% densities. In the next example, we show how one can use an equation of
% state (Peng-Robinson) to determine the equilibrium state instead. We
% start by constructing the necessary input parameters to the PR EoS for a
% system consisting of three named components. These properties are
% generated from CoolProp and we display them using a custom-made disp
% function.
eosname='soreide-whitson'; %SW model
%eosname='Peng-Robinson';
eosmodel = EquationOfStateModel([], compFluid, eosname);


% initialisation
z=[0.8 , 0.2];
patm=10^5;
CasTest=3;
switch CasTest
    case 1
        eosmodel.msalt=0;
        pres=[37.108,79.366,121.706,29.272,60.213,93.426]*barsa;
        Temp=[323.18,323.18,323.19,372.71,372.73,372.72]*Kelvin;
        xliqCO2Exp=[0.000461,0.001030,0.001544,0.000396,0.000857,0.001368];
    case 2
        eosmodel.msalt=0;
        Temp=[298.20,298.05,298.15,298.15,323.55,323.50,323.85,323.55,...
            323.30,323.35,373.85,373.80,373.85,373.65]*Kelvin;
        pres=[100.01, 150.01, 200.01, 200.01, 101.11, 101.31, 130.01,...
            165.01, 199.91, 200.11, 100.01, 100.01, 100.01, 175.11]*barsa;     
        xliqH2Exp=[0.00135994,0.00199679,0.00264397,0.00263320,0.00125091,...
            0.00123925,0.00159253,0.00201117,0.00244186,0.00245479,...
            0.00142471,0.00140332,0.00140893,0.00243347];

    case 3
        eosmodel.msalt=1;
        Temp=[298.20,298.30,298.15,298.30,323.20,323.40,323.35,323.20,...
            323.30,323.30,323.20,373.25,373.40,373.10,373.15,373.45,373.15,373.00];
        pres=[100.71,150.01,150.01,200.01,100.31,100.61,101.01,150.06,...
            175.01,199.91,200.01,100.11,126.01,150.36,150.46,150.71,175.51,200.46];
        xliqH2Exp=[0.00107012,0.00159298,0.00161244,0.00213575,0.00102099,...
            0.00102078,0.00102426,0.00151377,0.00176728,0.00202590,...
            0.00204487,0.00119671,0.00148492,0.00175595,0.00179573,...
            0.00177350,0.00204000,0.00234604];
    case 4
        eosmodel.msalt=2;
        Temp=[298.15,298.05,298.05,323.20,323.40,323.35,323.40,373.05,373.20,373.40];
        pres=[100.01,150.01,200.01,100.01,150.01,150.01,200.01,100.01,150.01,200.51];
        xliqH2Exp=[0.00088640,0.00132235,0.00171848,0.00088260,0.00131242,...
            0.00128912,0.00172402, 0.00099379, 0.00151866, 0.00205031 ];
    case 5
        eosmodel.msalt=4;
        Temp=[298.20, 298.20, 298.15, 323.30, 323.40, 323.40, 373.25, 373.35, 373.15];
        pres=[100.01, 150.01, 200.01, 100.01, 150.01, 200.01, 100.01, 150.01, 200.01];
        xliqH2Exp=[0.00059422, 0.00093595, 0.00121838, 0.00061736, ...
            0.00095752, 0.00129237, 0.00077991, 0.00114509, 0.00157469];
end

        

% Perform the standalone flash calculation to determine liquid fraction 
namecp=eosmodel.CompositionalMixture.names; %components name
nc=size(pres,2);
indH2=find(strcmp(namecp,'H2'));
indCO2=find(strcmp(namecp,'CO2'));
indH2O= find(strcmp(namecp,'H2O'));


[L, x, y, reports] = standaloneFlash(pres, Temp, z, eosmodel); 
if ~isempty(indH2)
    xliqH2=x(:,indH2);
    namefile=['H2solubility_cas',num2str(CasTest),'_msalt',num2str(eosmodel.msalt),eosname,'.dat'];
    fiD = fopen(namefile,'wt'); 
    for i=1:nc
        fprintf(fiD,'%6.2f %12.4f %12.8f %12.8f\n',Temp(i),pres(i)/patm,xliqH2(i),xliqH2Exp(i));
    end
    fclose(fiD);
end
if ~isempty(indCO2)
    xliqCO2=x(:,indCO2);
     namefile=['CO2solubility_cas',num2str(CasTest),'_msalt',num2str(eosmodel.msalt),eosname,'.dat'];
    fiD = fopen(namefile,'wt'); 
    for i=1:nc
        fprintf(fiD,'%6.2f %12.4f %12.8f %12.8f\n',Temp(i),pres(i)/patm,xliqCO2(i),xliqCO2Exp(i));
    end
    fclose(fiD);
end



% %write data in a file
% namefile=['H2solubility_cas',num2str(CasTest),'_msalt',num2str(eosmodel.msalt),eosname,'.dat'];
% fiD = fopen(namefile,'wt'); 
% %fprintf(fiD,'Temperature (K), Pressure (bar), H2 molar fraction \n');
% for i=1:nc
%    fprintf(fiD,'%6.2f %12.4f %12.8f %12.8f\n',Temp(i),pres(i)/patm,xliqH2(i),xliqH2Exp(i));
% end
% fclose(fiD);
