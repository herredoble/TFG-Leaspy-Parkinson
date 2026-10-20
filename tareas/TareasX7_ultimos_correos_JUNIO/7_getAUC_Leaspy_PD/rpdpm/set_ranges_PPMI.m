function [data_rpdpm,ranges]=set_ranges_PPMI(data_rpdpm)

% "MOCA"
% "GDS"
% "REM"
% "LNS"
% "SDMT" (5)
% "JLO"
% "INMEDIATE"
% "HVLTRDLY"
% "SFT"
% "SCOPA" (10)
% "MEANCAU"
% "ABETA"
% "TAU"
% "UPDRS3"

data_rpdpm.MOCA(data_rpdpm.MOCA<10)=10;
data_rpdpm.SDMT(data_rpdpm.SDMT>80)=80;
data_rpdpm.SFT(data_rpdpm.SFT>100)=100;
data_rpdpm.ABETA(data_rpdpm.ABETA>2000)=2000;
data_rpdpm.TAU(data_rpdpm.TAU>400)=400;
data_rpdpm.TAU(data_rpdpm.TAU<80)=80;
data_rpdpm.UPDRS3(data_rpdpm.UPDRS3>80)=80;

% Ranges of biomarkers
ranges = [
    10 30;... MOCA
    0 15;... GDS
    0 14;... REM
    0 20;... LNS
    0 80;... SDMT (5)
    0 15;... JLO
    0 40;... IMMEDIATE
    0 12;... HVLTRDLY
    7 100;... SFT
    0 45;... SCOPA (10)
    0.1650 3.7050;... MEANCAU
    200 2000;... ABETA
    80 400;...TAU
    0 80 ....UPDRS3
    ];

end