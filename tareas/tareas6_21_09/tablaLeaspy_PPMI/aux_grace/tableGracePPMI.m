function [grace,percent,point_feat,change_sign_feat]=tableGracePPMI(show,n_dataset)
clc;
if(show) 
%     close all
end
addpath('./aux_grace');
addpath('./external/lme/univariate');

population={'FS_507','PA_633','PPMISep20'};
dataset=population{n_dataset};

if(strcmp(dataset,'PPMISep20'))
    [PPMI_PD,PPMI_NC]=cognitivePD_PPMI;
    PPMI_PD.Month=PPMI_PD.YEAR*12;
    PPMI_NC.Month=PPMI_NC.YEAR*12;
    %Only for showing
    PPMI_NC.onset=(PPMI_NC.Month-(15*12))*30;
    PPMI_PD.onset=(PPMI_PD.Month+(15*12))*30;
    
    PPMI=[PPMI_NC;PPMI_PD];
    PPMI.fsidbase=PPMI.PATNO;
    %PPMI.Converts=PPMI.MCI;
    PPMI.Converts=PPMI.Convert;
    %PPMI.diagnose=categorical(PPMI.MCI,[0 1],{'Control','PD'});
    PPMI.diagnose=categorical(PPMI.Convert,[0 1],{'Control','PD'});
else
    [PPMI_NC,PPMI_PD]=selectedSubjectPPMI(dataset);
    %Only for showing
    PPMI_NC.onset=(PPMI_NC.Month-(15*12))*30;
    PPMI_PD.onset=(PPMI_PD.Month+(15*12))*30;
    PPMI=[PPMI_NC;PPMI_PD];

end

%Only for showing
group=[zeros(size(PPMI_NC,1),1);ones(size(PPMI_PD,1),1)];


%% Mean zero for study time
RID=unique(PPMI.fsidbase);
ID=zeros(size(PPMI,1),1);
t=ID;
for i=1:length(RID)
    index=find(PPMI.fsidbase==RID(i));
    ID(index)=i;
    months=PPMI.Month(index);
    time_zero_study=max(months)/2;
    t(index)=PPMI.Month(index)-time_zero_study;

end
PPMI.ID=ID;
PPMI.t=t;

%% Features
    
if(strcmp(dataset,'FS_507'))
    %RAV=PPMI.RAmigdala./PPMI.ICV*1e3;
    NCV=(PPMI.LCaudete+PPMI.RCaudete)./PPMI.ICV*1e3;
    SPT_RH=PPMI.R_SupPariet;

    %raw=[RAV,PPMI.UPDRSTOT,PPMI.GDS,PPMI.SDMTOTAL,PPMI.VLTANIM];
    %feature_names={'RAV','UPDRSTOT','GDS','SDMTOTAL','VLTANIM'};
    raw=[NCV,SPT_RH,PPMI.UPDRSTOT,PPMI.SDMTOTAL];
    feature_names={'NCV','SPTRH','UPDRSTOT','SDMTOTAL'};
    
%     numMRImarkers=2;
%     bMRImarker(1:numMRImarkers)=false;
    

elseif(strcmp(dataset,'PA_633'))
    %LAV=PPMI.LAmygdala./PPMI.ICV*1e3;
    %LPaV=PPMI.LPallidum./PPMI.ICV*1e3;
    IPT=PPMI.L_InfPariet+PPMI.R_InfPariet;
    %NCV=(PPMI.LCaudete+PPMI.RCaudete)./PPMI.ICV*1e3;
    %RCV=PPMI.RCaudete./PPMI.ICV*1e3;
    ITT_RH=PPMI.R_InfTemp;
    %SPT_RH=PPMI.R_SupPariet;
    ECT_RH=PPMI.R_ECT;
    %TPT=PPMI.L_TempPole+PPMI.R_TempPole;
    %RHV=PPMI.RHippVol./PPMI.ICV*1e3;
    %NHV=(PPMI.RHippVol+PPMI.LHippVol)./PPMI.ICV*1e3;

    raw=[IPT,ECT_RH,ITT_RH];
    feature_names={'IPT','ECT_RH','ITT_RH'};
    
    %raw=[LAV,IPT,PPMI.R_ECT,PPMI.R_InfTemp,LPaV];
    %feature_names={'LAV','IPT','ECT_RH','ITT_RH','LPaV'};
    %raw=[NCV,IPT,ECT_RH,ITT_RH];
    %feature_names={'NCV','IPT','ECT_RH','ITT_RH'};
    %raw=[RCV,SPT_RH,NHV];
    %feature_names={'RCV','SPT_RH','NHV'};
    
%     numMRImarkers=3;
%     bMRImarker(1:numMRImarkers)=false;
    

    
elseif(strcmp(dataset,'PPMISep20'))
    feature_names={'LNS','SDMT','HVLTRDLY','SFT','CONTCAU','ABetaTAU','UPDRS1'};
    raw=[PPMI.lns,PPMI.SDMTOTAL,PPMI.HVLTRDLY,PPMI.sft,PPMI.con_caudate,...
        PPMI.abeta./PPMI.tau,PPMI.updrs1_score];

%     feature_names={'LNS','SDMT','SFT','SCOPA','CONTCAU','ABETA'};
%     raw=[PPMI.lns,PPMI.SDMTOTAL,PPMI.sft,PPMI.scopa,PPMI.con_caudate,...
%         PPMI.abeta];
    
    change_sign_feat=true(size(raw,2),1);
    change_sign_feat(end)=false; %UPDRS1
    PPMI.years=PPMI.YEAR;
    
    
else
end

%% Marker correction
corrAge=true;
if(corrAge)
    numCorrMarkers=size(raw,2);
    bMRImarker(1:numCorrMarkers)=false;
    ni_all=get_ni(PPMI.years,length(unique(PPMI.fsidbase)));
    SexEduGrp=false; %only age with MRI markers
    bTwoGroup=false; %for correction without clinical groups
    index_NC=1:size(PPMI_NC,1);

    LMEmodels=getLME_univariate_PPMI(PPMI(index_NC,:),raw(index_NC,1:numCorrMarkers),...
        feature_names(1:numCorrMarkers),bMRImarker,SexEduGrp,bTwoGroup);
    % features_corr=correctionCovariablesLME_PPMI(LMEmodels,PPMI,raw(:,1:numCorrMarkers),ni_all);
    features=correctionCovariablesLME_PPMI(LMEmodels,PPMI,raw(:,1:numCorrMarkers),ni_all);
else
    features=raw;
end

% if(strcmp(dataset,'FS_507'))
%     features=[-features_corr,raw(:,3),-raw(:,4)];
% elseif(strcmp(dataset,'PA_633'))
%     %features=[-features_corr(:,1:2),features_corr(:,3)]; %IPT,R_ECT,R_InfTemp
%     %features=[-features_corr(:,1:3),features_corr(:,4)]; %NVC,IPT,R_ECT,R_InfTemp
%     features=-features_corr; 
% else
% end

%% Incresead features
for i=1:length(change_sign_feat)
    if(change_sign_feat(i)==true)
        features(:,i)=-features(:,i);
    end
end
%% Transforming
%[features,~,~]=percentil_interp(features,group);
[features,percent,point_feat]=percentil_interp(features,group);

if(show)
    level_linear=0;

    
    for i=1:size(features,2)
        figure(i);clf;hold on;
        title(feature_names{i});
        mask_nan=~isnan(features(:,i));

        draw_longTrayPPMI(PPMI.fsidbase(mask_nan),PPMI.onset(mask_nan),...
             features(mask_nan,i),PPMI.diagnose(mask_nan),level_linear);
    end
    
end



%% Table for R
num_features=size(features,2);

grace=zeros(length(features(:)),7);

index_begin=1;
for i=1:length(RID)
    index=find(PPMI.ID==i);
    index_end=index_begin+(length(index)*num_features)-1;
    features_subj=features(index,:);
    raw_subj=raw(index,:);
    id_marker=zeros(size(features_subj));
    for j=1:num_features
        id_marker(:,j)=j;
    end
    [t,index_time]=sort(PPMI.t(index));
    features_subj=features_subj(index_time,:);
    raw_subj=raw_subj(index_time,:);
    grace(index_begin:index_end,1)=repmat(t/12,num_features,1);    
    grace(index_begin:index_end,2)=features_subj(:);
    grace(index_begin:index_end,3)=raw_subj(:);   
    grace(index_begin:index_end,4)=id_marker(:);
    grace(index_begin:index_end,5)=repmat(PPMI.ID(index),num_features,1);
    grace(index_begin:index_end,6)=repmat(PPMI.Converts(index),num_features,1);
    grace(index_begin:index_end,7)=repmat(PPMI.fsidbase(index),num_features,1);
    grace(index_begin:index_end,8)=repmat(PPMI.convertTime(index),num_features,1);
    
    index_begin=index_end+1;
    
end

grace=grace(~isnan(grace(:,2)),:);
grace=array2table(grace,'VariableNames',...
    {'t','Y','raw','Outcome','id','Convert','RID','ConvertTime'});

% writetable(grace,'dd.xlsx');

%% Correlations
if(show)
    outcomes=unique(grace.Outcome);
    LMEmodels=cell(length(outcomes),1);
    grace=cleanID_sinVisit(grace);

    for i=1:length(outcomes)
        grace_outcome=grace(grace.Outcome==outcomes(i),:);
        LMEmodels{i}=getLME_outcome(grace_outcome);
    end
    [inter,slope]=corrOutcomeLME(LMEmodels);

    % corr_inter=corrcov(cov(inter));
    % corr_slope=corrcov(cov(slope));
    
    figure(size(features,2)+1);
    corrplot(array2table(inter,'VariableNames',feature_names));
    figure(size(features,2)+2);
    corrplot(array2table(-slope,'VariableNames',feature_names));
end    

%% exit
% rmpath('./aux_grace');
rmpath('./external/lme/univariate');


end



