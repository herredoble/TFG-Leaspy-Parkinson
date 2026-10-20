function [data_rpdpm,ranges]=set_ranges_PPMI_LSHSM(data_rpdpm)

data_rpdpm.SDMT(data_rpdpm.SDMT>80)=80;
data_rpdpm.ABETA(data_rpdpm.ABETA>2000)=2000;
data_rpdpm.SFT(data_rpdpm.SFT>100)=100;

% Ranges of biomarkers
%  ABETA, HVLTRDLY, LNS, MEANCAU, SDMT,SFT
ranges = [200 2000;...ABETA
    0 12;...HVLTRDLY
    0 20;...LNC
    0.1650 3.7050;... MEANCAU
    0 80;... SDMY
    7 100;... SFT
    ];

end