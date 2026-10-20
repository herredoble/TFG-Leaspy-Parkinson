function [PPMI,features,features_names,ni,features_SNan] = getDataPPMISep20(n_feat,n_cohort)
% This function is charged of importing and cleaning the Tresults
% table, taking only candidate features. However, Tresults is equally
% returned so cleaned variables are still accessible. Extracted
% features can be Rois, neuropsycological tests, clusters
% or all of them.
%

%% Data input
cohort={'all','onlyPD','Abeta_PD','onlyHC'};
typeFeatures={'NM','DAT','CSF','MOTOR','all'};    

if(strcmp(cohort{n_cohort},'all'))
    load('./data/PPMIApr20_5Years','PPMIApr20');
    mask_PD=PPMIApr20.APPRDX==1;
    mask_HC=PPMIApr20.APPRDX==2;
    
    PPMI=PPMIApr20(mask_PD |mask_HC,:);
    clear PPMIApr20
%     mask_PD=PPMI.APPRDX==1;
%     mask_HC=PPMI.APPRDX==2;


%% PD_MCI vs PD_HC
elseif(strcmp(cohort{n_cohort},'onlyPD'))
    [PD_MCI,PD_HC]=cognitivePD_PPMI;
    PPMI=[PD_MCI;PD_HC];
%     mask_PD=false(size(PPMI,1),1);
%     mask_HC=mask_PD;
%     mask_PD(1:size(PD_MCI,1))=true;
%     mask_HC(size(PD_MCI,1)+1:end)=true;
%% PD_AB+ vs PD_AB-
elseif(strcmp(cohort{n_cohort},'Abeta_PD'))
    [PD_ABmas,PD_ABmenos]=Abeta_PD_PPMI;
    PPMI=[PD_ABmas;PD_ABmenos];
%     mask_PD=false(size(PPMI,1),1);
%     mask_HC=mask_PD;
%     mask_PD(1:size(PD_ABmas,1))=true;
%     mask_HC(size(PD_MPD_ABmasCI,1)+1:end)=true;    
else
end    
    
   
%% BASICS
% These features are useful in the function.
time_visit=PPMI.YEAR;
ni = get_ni(PPMI.YEAR,numel(unique(PPMI.PATNO)));

%% Cognitive markers
NM=[PPMI.moca,PPMI.gds,PPMI.upsit,PPMI.rem,...
    PPMI.lns,PPMI.SDMTOTAL,PPMI.bjlot,PPMI.hvlt_immediaterecall,PPMI.HVLTRDLY,...
    PPMI.sft,PPMI.MSEADLG,PPMI.scopa];
NM_SNan = clean_invalid_CAM(NM,ni,time_visit);
NM_names={'MOCA','GDS','UPSIT','REM','LNS','SDMT','JLO','INMEDIATE','HVLTRDLY','SFT','ADL','SCOPA'};

%% DAT imaging
DatScan=[PPMI.con_caudate,PPMI.ips_caudate,PPMI.con_putamen,PPMI.ips_putamen,...
    PPMI.mean_putamen,PPMI.mean_caudate,PPMI.ai_putamen,PPMI.ai_caudate];
DatScan_SNan = clean_invalid_CAM(DatScan,ni,time_visit);
DatScan_name={'CONTCAU','IPSCAU','CONPUT','IPSPUT','MEANPUT','MEANCAU','AIPUT','AICAU'};

%% CSF
CSF=[PPMI.abeta,PPMI.asyn,PPMI.tau,PPMI.ptau,PPMI.abeta./PPMI.tau];
CSF_SNan = clean_invalid_CAM(CSF,ni,time_visit);
CSF_name={'ABETA','ASYN','TAU','PTAU','ABetaTAU'};

%% MOTOR
MOTOR=[PPMI.updrs_totscore,PPMI.updrs_totscore_on,PPMI.updrs1_score,...
    PPMI.updrs2_score,PPMI.updrs3_score,PPMI.updrs3_score_on];
MOTOR_SNan = clean_invalid_CAM(MOTOR,ni,time_visit);
MOTOR_name={'UPDRSTOT','UPDRSTOTON','UPDRS1','UPDRS2','UPDRS3','UPDRS3ON'};

%% Appropriate combination
if strcmp(typeFeatures{n_feat},'NM')
    features = NM;
    features_SNan=NM_SNan;
    features_names = NM_names;
elseif strcmp(typeFeatures{n_feat},'DAT')
    features = DatScan;
    features_SNan=DatScan_SNan;
    features_names = DatScan_name;
elseif strcmp(typeFeatures{n_feat},'CSF')
    features = CSF;
    features_SNan=CSF_SNan;
    features_names = CSF_name;
elseif strcmp(typeFeatures{n_feat},'MOTOR')
    features = MOTOR;
    features_SNan=MOTOR_SNan;
    features_names = MOTOR_name;

elseif strcmp(typeFeatures{n_feat},'all')
    features = [NM,DatScan,CSF,MOTOR];
    features_SNan=[NM_SNan,DatScan_SNan,CSF_SNan,MOTOR_SNan];
    features_names = [NM_names,DatScan_name,CSF_name,MOTOR_name];

else
    error('Invalid input statement. feat_sel must be ROI, CAM, MASS or all');
end

end




