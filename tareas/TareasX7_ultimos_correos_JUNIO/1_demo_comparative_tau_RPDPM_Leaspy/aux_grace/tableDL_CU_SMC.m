function tbl_DL=tableDL_CU_SMC(n_typeFeat,n_trial)

close all
clc;


%addpath('./aux_grace');
%addpath('./external/lme/univariate');


typeFeatures={'ROI','CAM','ROICAM','all'};
type_features=cellstr(typeFeatures{n_typeFeat});
%% ADNI data
[ADNI,features,feature_names,~,features_SNan] = getDataADNIMERGE_v1_2b(type_features);

% mask_CN=ADNI.stMCI==0; % CU
mask_CN=(ADNI.stMCI==0) | (ADNI.stMCI<0); % CU and SMC
% mask_sMCI=ADNI.stMCI==1 & ADNI.Convert==0;
% mask=mask_CN | mask_sMCI;

mask= mask_CN & ADNI.M==0;

ADNI=ADNI(mask,:);
features=features(mask,:);
features_SNan=features_SNan(mask,:);
feature_names=string(feature_names);

group=ADNI.Convert;
% group(group<0)=0; %SMC

%% Amyloid positive
% th_pTAU_AB=0.028;
% pTAU_AB=ADNI.PTAU./ADNI.ABETA;
% ADNI.Apos= zeros(size(ADNI,1),1);
% mask_no_nan=~isnan(pTAU_AB);
% %Apos: 0 (nan) 1 (A-) 2 (A+)
% ADNI.Apos(mask_no_nan)=double(pTAU_AB(mask_no_nan)>=th_pTAU_AB)+1;





%% Mean zero for study time
% RID=unique(ADNI.RID);
% t=zeros(size(ADNI,1),1);
% for i=1:length(RID)
%     index=find(ADNI.RID==RID(i));
%     months=ADNI.M(index);
%     time_zero_study=max(months)/2;
%     t(index)=ADNI.Month(index)-time_zero_study;
% 
% end
% ADNI.ID=ADNI.RID;
% ADNI.t=t;

%% Features
% Marcadores CAM:
    % 01 - RAVLT_forgetting (+)
    % 02 - RAVLT_immediate (-) 
    % 03 - RAVLT_learning (-)
    % 04 - RAVLT_perc_forgetting(+)
    % 05 - ADAS13 (+)
    % 06 - FAQ (+)
    % 07 - MMSE (-)
    % 08 - CDRSB (+)
    % 09 - LDELTOTAL (-) 
    % 10 - PACCdigit (-)
    % 11 - PACCtrailsB (-) 
    % 12 - TRABSCOR (+)


    
if(strcmp(type_features,'CAM'))
    if(n_trial==1)
        n_CorrAge=3; %SinAge
        selectedFeatures=9;%LDELTOTAL
        %numMRImarkers=0;
        change_sign_feat=true(length(selectedFeatures),1);
    elseif(n_trial==2)
        n_CorrAge=3; %SinAge
        selectedFeatures=[3,5,8,12];% Learning- A13+ CDRSB+ TRABSCOR+
        %numMRImarkers=2;
        change_sign_feat=false(length(selectedFeatures),1);
        change_sign_feat(1)=true; %Learning-                    
    elseif(n_trial==3)
        n_CorrAge=2; %Age
        selectedFeatures=[3,5,8,12];% Learning- A13+ CDRSB+ TRABSCOR+
        change_sign_feat=false(length(selectedFeatures),1);
        change_sign_feat(1)=true; %Learning-                    
    elseif(n_trial==4)
        n_CorrAge=2; %Age
        selectedFeatures=[5,8];% A13+ CDRSB+
        change_sign_feat=false(length(selectedFeatures),1);
    elseif(n_trial==5)
        n_CorrAge=3; %SinAge
        selectedFeatures=[5,8,12];% A13+ CDRSB+ TRABSCOR+
        change_sign_feat=false(length(selectedFeatures),1);
    elseif(n_trial==6)
        n_CorrAge=3; %SinAge
        selectedFeatures=[3,5,6,8,11];% Learning- A13+ FAQ+ CDRSB+ PACCTrial-
        %numMRImarkers=2;
        change_sign_feat=false(length(selectedFeatures),1);
        change_sign_feat([1,5])=true; %Learning- PACCTrial-  
    elseif(n_trial==7)
        n_CorrAge=2; %Age
        selectedFeatures=[2,3,5,8,12];% Immediate- Learning- A13+ CDRSB+ TRABSCOR+
        change_sign_feat=false(length(selectedFeatures),1);
        change_sign_feat([1,2])=true; %Immediate- Learning-                    
        
    else
    end

else
end
    
    

% raw=features(:,selectedFeatures);
raw=features_SNan(:,selectedFeatures);

feature_names=feature_names(selectedFeatures);
% years=ADNI.M/12;
ADNI.years=ADNI.M/12;

% bMRImarker=false(length(selectedFeatures),1);





%% Marker correction
typeCov={'CorrAge','AddAge','SinAge','AddAgeEduc'};
corrAge=typeCov{n_CorrAge};

if(strcmp(corrAge,'CorrAge'))
%     ni_all=get_ni(years,length(unique(ADNI.RID)));
%     fld=getFLD(ADNI.FLDSTRENG,ADNI.FLDSTRENG_bl,ADNI.ORIGPROT,...
%         ni_all,ADNI.RID)-1;
%     ADNI.fld=fld;
% 
%     SexEduGrp=false; %only age and fld with MRI markers
%     bTwoGroup=false; %for correction without clinical groups
%     mask_sCN= (ADNI.stMCI == 0) & (ADNI.Convert == 0); 
%     LMEmodels=getLME_univariate_ADNI(ADNI(mask_sCN,:),raw(mask_sCN,:),...
%         feature_names,bMRImarker,SexEduGrp,bTwoGroup);
%     features=correctionCovariablesLME_ADNI(LMEmodels,ADNI,raw,ni_all);
elseif(strcmp(corrAge,'AddAge'))
   raw=[raw,ADNI.years+ADNI.AGE];
   features=raw;
   change_sign_feat=[change_sign_feat;false];
   feature_names=[feature_names;'Age'];
elseif(strcmp(corrAge,'AddAgeEduc'))
%    raw=[raw,ADNI.years+ADNI.AGE,ADNI.PTEDUCAT];
%    features=raw;
%    change_sign_feat=[change_sign_feat;false;true];
%    feature_names=[feature_names;'Age';'Education'];   
   
else
    features=raw;
end


%% Incresead features
% for i=1:length(change_sign_feat)
%     if(change_sign_feat(i)==true)
%         features(:,i)=-features(:,i);
%     end
% end
%% Transforming
% [features,percent,point_feat]=percentil_interp(features,group);


%% Table for R
colnames=[feature_names;"group"];
tbl_DL=array2table([features,group],'VariableNames',colnames);

tbl_DL.group=categorical(tbl_DL.group,[0 1],{'stable','convert'});


%% exit
rmpath('./aux_grace');
% rmpath('./external/lme/univariate');


end



