function [PD_mas,PD_menos]=Abeta_PD_PPMI

%addpath('./selected_aux');

load('./data/PPMIApr20_5Years','PPMIApr20');

mask_PD=PPMIApr20.APPRDX==1;
PPMI_PD=PPMIApr20(mask_PD,:);

% PPMI_PD=CSF_groups(PPMI_PD);
% PD_mas=PPMI_PD(PPMI_PD.TAUABETAmas==1,:);
% PD_menos=PPMI_PD(PPMI_PD.TAUABETAmas==0,:);


%% AB_bsl
PPMI_PD=PPMI_PD(PPMI_PD.YEAR==0,:);
PPMI_PD=PPMI_PD(~isnan(PPMI_PD.abeta),:);
thres_Abeta=880; % Hasson18
mask_ABmas_bsl=PPMI_PD.abeta<thres_Abeta;
PD_mas=PPMI_PD(mask_ABmas_bsl==1,:);
PD_menos=PPMI_PD(mask_ABmas_bsl==0,:);


%% delete
%rmpath('./selected_aux');

end





