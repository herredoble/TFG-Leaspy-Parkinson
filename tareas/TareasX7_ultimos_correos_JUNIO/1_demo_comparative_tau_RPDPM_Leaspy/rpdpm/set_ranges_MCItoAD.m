function [data_rpdpm,ranges]=set_ranges_MCItoAD(data_rpdpm)

data_rpdpm.FORGETTING(data_rpdpm.FORGETTING<-5)=-5;
data_rpdpm.FORGETTING(data_rpdpm.FORGETTING>15)=15;
%data_rpdpm.IMMEDIATE (data_rpdpm.IMMEDIATE<5)=5;
data_rpdpm.PERCFORGETTING (data_rpdpm.PERCFORGETTING<0)=0;
%data_rpdpm.FAQ (data_rpdpm.FAQ>15)=15;
%data_rpdpm.MMSE (data_rpdpm.MMSE<15)=15;
data_rpdpm.PACCTRAILSB (data_rpdpm.PACCTRAILSB<-40)=-40;
data_rpdpm.TRABSCOR (data_rpdpm.TRABSCOR>300)=300;


% Ranges of biomarkers
ranges = [...
    -5 15;... FORGETTING
    0 80;... IMMEDIATE
    -8 14;... LEARNING
    0 100;... PERCFORGETTING
    0 85;...ADAS13
    0 30;...FAQ
    0 30;...MMSE
    0 18;... CDRSB
    0 30;... LDELTOTAL
    -40 10;...PACCTrial
    0 300 %TRABSCOR
    ];

end