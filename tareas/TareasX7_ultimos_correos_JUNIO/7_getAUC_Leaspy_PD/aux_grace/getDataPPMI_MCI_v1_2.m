function [PPMI,features,features_names,ni,features_SNan] = getDataPPMI_MCI_v1_2(n_feat,n_cohort)
% This function is charged of importing and cleaning the Tresults
% table, taking only candidate features. However, Tresults is equally
% returned so cleaned variables are still accessible. Extracted
% features can be Rois, neuropsycological tests, clusters
% or all of them.
%

%% Data input
cohort={'PPMIApr20_5Years'};
typeFeatures={'NM','DAT','CSF','MOTOR','all'};    


%% PD_MCI vs PD_HC
if(strcmp(cohort{n_cohort},'PPMIApr20_5Years'))
    load('./data/PPMIApr20_5Years','PPMIApr20');
    mask_PD=PPMIApr20.APPRDX==1;
    PPMI=PPMIApr20(mask_PD,:);
    PPMI=sortrows(PPMI,{'PATNO','YEAR'});
    
    [cogstate_interp,ID_sinCog] = fill_cogstate_with_NAN(PPMI);
    PPMI.cogstate=cogstate_interp;
    [PPMI,ID_sinCog]=convertTime_PPMI(PPMI);
    %PPMI=estimateCogState(PPMI,ID_sinCog);
    PPMI.date=convert_visit_date_PPMI(PPMI.visit_date);
    PPMI.date_bsl=getDate_bsl(PPMI.PATNO,PPMI.date,PPMI.YEAR);
    PPMI.onset=getOnsetTime_PPMI(PPMI.PATNO,PPMI.date,...
        PPMI.cogstate_end>1,PPMI.convertTimeCogstate,PPMI.YEAR);
    
 
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
% DatScan=[PPMI.con_caudate,PPMI.ips_caudate,PPMI.con_putamen,PPMI.ips_putamen,...
%     PPMI.mean_putamen,PPMI.mean_caudate,PPMI.ai_putamen,PPMI.ai_caudate];
DatScan=[PPMI.mean_putamen,PPMI.mean_caudate];

DatScan_SNan = clean_invalid_CAM(DatScan,ni,time_visit);
% DatScan_name={'CONTCAU','IPSCAU','CONPUT','IPSPUT','MEANPUT','MEANCAU','AIPUT','AICAU'};
DatScan_name={'MEANPUT','MEANCAU'};


%% CSF
CSF=[PPMI.abeta,PPMI.asyn,PPMI.tau,PPMI.ptau,PPMI.ptau_ab];
CSF_SNan = clean_invalid_CAM(CSF,ni,time_visit);
CSF_name={'ABETA','ASYN','TAU','PTAU','pTauAB'};

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




