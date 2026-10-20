function [Leaspy_train,Leaspy_test]=tableLeaspy_PPMI

clc;
n_dataset=1;
n_cog_criteria=3;
n_trial=0; % SDMT-
% n_trial=18; % LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-

addpath('./aux_grace');

population={'PPMISep20'};
dataset=population{n_dataset};

cognitiveCriteria={'MCI_test','UPDRS1','SITE'};
varCogCriteria=cellstr(cognitiveCriteria{n_cog_criteria});


if(strcmp(dataset,'PPMISep20'))
    n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
    %[PPMI,features,feature_names,~,features_SNan] = getDataPPMI_MCI_v1(n_typeFeat,n_dataset);
    [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_MCI_v1_1(n_typeFeat,n_dataset);
    
    % Cognitive criteria
    if(strcmp(varCogCriteria,'SITE'))
        % PPMI.Convert=double(PPMI.cogstate_end>1);
        PPMI.Convert=PPMI.cogstate_end;
        PPMI.convertTime=PPMI.convertTimeCogstate;
    elseif(strcmp(varCogCriteria,'MCI_test'))
        PPMI.Convert=double(PPMI.MCI_testscores_end);
        PPMI.convertTime=PPMI.convertTimeMCI_test;
    elseif(strcmp(varCogCriteria,'UPDRS1'))
        PPMI.Convert=double(PPMI.np1cog_end>1);
        PPMI.convertTime=PPMI.convertTimeNp1cog;
    else
    end
   
    PPMI.Month=PPMI.YEAR*12;
    feature_names=string(feature_names)';
end

%% Features
mask_bsl=PPMI.YEAR==0;
feat_bsl=features(mask_bsl,:);
feat_bsl_SNan=features_SNan(mask_bsl,:);
mask_bsl_nan=isnan(feat_bsl);
feat_bsl(mask_bsl_nan)=feat_bsl_SNan(mask_bsl_nan);
features(mask_bsl,:)=feat_bsl;

if(strcmp(dataset,'PPMISep20'))
    
    if(strcmp(cognitiveCriteria{n_cog_criteria},'SITE'))
        if(n_trial==0)
            n_CorrAge=3;%SinAge
            selectedFeatures=6;%SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1);       
        elseif(n_trial==1) %ProposedAbril21
            n_CorrAge=2;%AddAge
            %selectedFeatures=[4,5,6,9,10,18,25];%REM LNS SDMT HVLTRDLY SFT MEANCAU ABetaTAU (CogState2)
            selectedFeatures=[4,5,6,9,10,14,19];%REM LNS SDMT HVLTRDLY SFT MEANCAU ABetaTAU (CogState2)
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(1)=false;%REM
        elseif(n_trial==2)
            n_CorrAge=3;%SinAge
            selectedFeatures=[5,6,10,12,18,25]; % LNS-, SDMT-,  SFT-, SCOPA+, CAU-, ABetaTAU-
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(4)=false;%SCOPA
         elseif(n_trial==3) %Schrag17
            n_CorrAge=2;%AddAge
            selectedFeatures=[3,4,14,15]; %UPSIT- REM+ MEANCAU- ABeta- (Schrag17)
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(2)=false;%REM
         elseif(n_trial==4)%Wilson20
            n_CorrAge=2;%AddAge
            selectedFeatures=[6,10]; % SDMT- SFT-
            change_sign_feat=true(length(selectedFeatures),1);
         elseif(n_trial==5) %AnovaDic21
            n_CorrAge=3;%SinAge
            selectedFeatures=[1,4:6,8,10,12,14,19]; % MOCA-, REM+, LNS-, SDMT-, IMM-, SFT-, SCOPA+, MEANCAU-,ABtTAU-
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat([2,7])=false;%REM, SCOPA
         elseif(n_trial==6) %Tang21
            n_CorrAge=3;%SinAge
            selectedFeatures=[1,4:6,8,9,20]; % MOCA-, REM+, LNS-, SDMT-, IMM-, DLY-, UPDRS+
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat([2,7])=false;%REM UPDRS
         elseif(n_trial==7)%Wilson20+Caspell_17
            n_CorrAge=2;%AddAge
            selectedFeatures=[6,10,14,15]; % SDMT- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1);
         elseif(n_trial==8)%Wilson20+Caspell_17
            n_CorrAge=2;%AddAge
            selectedFeatures=[4,5,6,10,14,15]; % REM+ LNS- SDMT- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1); 
            change_sign_feat(1)=false;%REM
         elseif(n_trial==9)%Wilson20+Caspell_17 (B)
            n_CorrAge=2;%AddAge
            selectedFeatures=[4,5,6,8,10,12,14,15]; 
            % REM+ LNS- SDMT- IMM- SFT- SCOPA+ MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1); 
            change_sign_feat([1,6])=false;%REM, SCOPA
         elseif(n_trial==10) %Dic21
            n_CorrAge=3;%SinAge
            selectedFeatures=[2,4,5,6,9,10,14,15,17];% GDS+ REM+ LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA- TAU+ 
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat([1,2,9])=false;%GDS, REM TAU
         elseif(n_trial==11) %Dic21
            n_CorrAge=3;%SinAge
            selectedFeatures=[2,4,5,6,9,10,14,19];% GDS+ REM+ LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA/TAU-
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat([1,2])=false;%GDS, REM            
         elseif(n_trial==12) %Tsiourus+Smith
            n_CorrAge=2;%AddAge
            selectedFeatures=[5,6,8,10,12]; % LNS- SDMT- IMM- SFT- SCOPA+(Tsiourus+Smith)
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(5)=false;%SCOPA+            
         elseif(n_trial==13)%ene22
            n_CorrAge=2;%AddAge
            selectedFeatures=[5,6,8,10]; %LNS- SDMT- IMM- SFT-
            change_sign_feat=true(length(selectedFeatures),1);
         elseif(n_trial==14) %Tsiourus+Smith
            n_CorrAge=2;%AddAge
            selectedFeatures=[5,6,9,12]; % LNS- SDMT- HVLTRDLY- SCOPA+ (Tsiourus+Smith B)
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(4)=false;%SCOPA+            
         elseif(n_trial==15) %Mar22
            n_CorrAge=2;%AddAge
            selectedFeatures=[5,6,7,8,12]; % LNS- SDMT- JLO- IMM- SCOPA+
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(5)=false;%SCOPA+             
          elseif(n_trial==16) %Mar22
            n_CorrAge=3;%SinAge
            selectedFeatures=[1,6,7,9];% MOCA- SDMT- JLO- HVLTRDLY-
            change_sign_feat=true(length(selectedFeatures),1); %PACC          
           elseif(n_trial==17) %Mar22
            n_CorrAge=2;%ConAge
            selectedFeatures=[1,4,8,9];% MOCA- REM+ IMM- HVLTRDLY-
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(2)=false;%REM+ 
           elseif(n_trial==18) %jpalomino
            n_CorrAge=3;%SinAge
            selectedFeatures=[5,6,9,10,14,15];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1);
           elseif(n_trial==19) %mayo22
            n_CorrAge=2;%AddAge
            selectedFeatures=[6,9,10,14,19];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- pTauAB+
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(5)=false;
           elseif(n_trial==20) %mayo22
            n_CorrAge=3;%SinAge
            selectedFeatures=[6,9,10,14,19];% SDMT- HVLTRDLY- SFT- MEANCAU- pTauAB+
            change_sign_feat=true(length(selectedFeatures),1);
            change_sign_feat(5)=false;
           elseif(n_trial==21) %sept22
            n_CorrAge=3;%SinAge
            selectedFeatures=[5,6,9,10];% LNS- SDMT- HVLTRDLY- SFT-
            change_sign_feat=true(length(selectedFeatures),1);
           elseif(n_trial==22) %sept22
            n_CorrAge=2;%AddAge
            selectedFeatures=[6,9,10];% SDMT- HVLTRDLY- SFT-
            change_sign_feat=true(length(selectedFeatures),1);
           elseif(n_trial==23) %see 18
            n_CorrAge=2;%SinAge
            selectedFeatures=[5,6,9,10,14,15];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1);
           elseif(n_trial==24) %see 20
            n_CorrAge=3;%SinAge
            selectedFeatures=[6,9,10,14,15];% SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
            change_sign_feat=true(length(selectedFeatures),1);
            
                        
        else
        end
    
        
    else
        
    end
    
    raw=features(:,selectedFeatures);    
    feature_names=feature_names(selectedFeatures);


    bMRImarker=false(length(selectedFeatures),1);
    PPMI.years=PPMI.YEAR;
    PPMI.gen=PPMI.gen-1;
    
else
end


    
    

%% Marker correction
typeCov={'CorrAge','AddAge','SinAge'};
corrAge=typeCov{n_CorrAge};

if(strcmp(corrAge,'CorrAge'))
%     SexEduGrp=false; %only age and fld with MRI markers
%     bTwoGroup=false; %for correction without clinical groups
%     %mask_stable= PPMI.Convert == 0; %???
%     mask_stable= PPMI.Convert == 1; %Cogstate
% 
%     LMEmodels=getLME_univariate_PPMI(PPMI(mask_stable,:),raw(mask_stable,:),...
%         feature_names,bMRImarker,SexEduGrp,bTwoGroup);
%     features=correctionCovariablesLME_PPMI(LMEmodels,PPMI,raw,ni);
    
elseif(strcmp(corrAge,'AddAge'))
   raw=[raw,PPMI.YEAR+PPMI.age];
   features=raw;
   change_sign_feat=[change_sign_feat;false];
   feature_names=[feature_names;'Age'];
   
else
    features=raw;
end




%% Incresead features
for i=1:length(change_sign_feat)
    if(change_sign_feat(i)==true)
        features(:,i)=-features(:,i);
    end
end
%% Transforming
[features,percent,point_feat]=percentil_interp(features,PPMI.Convert>1);


%% Table for Leaspy
TIME=PPMI.age+PPMI.YEAR;
Leaspy=[PPMI.PATNO,TIME,features,raw,PPMI.Convert,...
    PPMI.convertTime,PPMI.YEAR];
Leaspy=rmmissing(Leaspy);
colnames=["ID";"TIME";feature_names;feature_names+"_orig";"Convert";"ConverTime";"YEAR"];



% Cross varidation (train: 70%, test: 30%)
RID=unique(Leaspy(:,1));
cv = cvpartition(length(RID),'HoldOut',0.3);
idx = cv.test;
RID_train=RID(~idx);
RID_test=RID(idx);

mask_train=false(size(Leaspy,1),1);
for i=1:length(RID_train)
    index=find(Leaspy(:,1)==RID_train(i));
    if(length(index)>1)
         mask_train(index)=true;
     end
end
    
mask_test=false(size(Leaspy,1),1);
for i=1:length(RID_test)
    index=find(Leaspy(:,1)==RID_test(i));
    if(length(index)>1)
         mask_test(index)=true;
     end
end

Leaspy_train=Leaspy(mask_train,:);
Leaspy_test =Leaspy(mask_test,:);

Leaspy_train=rmmissing(Leaspy_train);
Leaspy_test=rmmissing(Leaspy_test);



Leaspy_train=array2table(Leaspy_train,'VariableNames',colnames);
Leaspy_test=array2table(Leaspy_test,'VariableNames',colnames);


if(n_trial==0)
    if(isfile('./data/SDMT_train.csv'))
         delete('./data/SDMT_train.csv');
    end
    writetable(Leaspy_train,'./data/SDMT_train.csv','Delimiter',',');
    
    if(isfile('./data/SDMT_test.csv'))
         delete('./data/SDMT_test.csv');
    end
    writetable(Leaspy_test,'./data/SDMT_test.csv','Delimiter',',');

elseif(n_trial==18)
    if(isfile('./data/LSDSMA_train.csv'))
         delete('./data/LSDSMA_train.csv');
    end
    writetable(Leaspy_train,'./data/LSDSMA_train.csv','Delimiter',',');
    
    if(isfile('./data/LSDSMA_test.csv'))
         delete('./data/LSDSMA_test.csv');
    end
    writetable(Leaspy_test,'./data/LSDSMA_test.csv','Delimiter',',');
end


%% exit
rmpath('./aux_grace');


end



