function summaryPPMI_PD_features_v1(n_Population)
n_Population=1
n_cog_criteria=3;
clc;
addpath('./selected_aux');

%% PD
cognitiveCriteria={'MCI_test','UPDRS1','SITE'};
varCogCriteria=cellstr(cognitiveCriteria{n_cog_criteria});

population={'PPMISep20_PD','PPMISep20_PD_HC'};
typeDataset=cellstr(population{n_Population});

n_typeFeat=5; %typeFeatures={'NM','DAT','CSF','MOTOR','all'};
n_cohort=1;

if(strcmp(typeDataset,'PPMISep20_PD'))
    [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_MCI_v1_1(n_typeFeat,n_cohort);
elseif(strcmp(typeDataset,'PPMISep20_PD_HC'))
    [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_HC_PD(n_typeFeat,...
        n_cohort);
else
end
    
% Cognitive criteria
if(strcmp(varCogCriteria,'SITE'))
    PPMI.MCI=PPMI.cogstate_end;
    PPMI.convertTime=PPMI.convertTimeCogstate;
    if(strcmp(typeDataset,'PPMISep20_PD'))
         group=categorical(PPMI.MCI,[1 2 3],{'PD-CU' 'PD-MCI' 'PD-Dementia'});
    elseif(strcmp(typeDataset,'PPMISep20_PD_HC'))
        group=categorical(PPMI.APPRDX,[1 2],{'PD' 'HC'});
    else
    end
else
end

    
        


if(n_typeFeat==5)
    fprintf('feature %s does not converge with LME. Only values are at baseline\n',string(feature_names(3)));
    %selectedFeatures=[1:2,4:size(features,2)]; % all
    selectedFeatures=1:size(features,2); % all
end

%% selectedFeatures
features=features(:,selectedFeatures);
features_SNan=features_SNan(:,selectedFeatures);
feature_names=string(feature_names);
feature_names=feature_names(selectedFeatures)';

% mask_bsl=PPMI.EVENT_ID=='BL';
mask_bsl=PPMI.YEAR==0';


%% Checking
features_SNan_bsl=features_SNan(mask_bsl,:);
for i=1:size(features_SNan_bsl,2)
    num_nan_feat=sum(isnan(features_SNan_bsl(:,i)));
    if(num_nan_feat>0)
        fprintf(2,'Subjects %d without %s marker\n',num_nan_feat,string(feature_names(i)));
    end
end

%% stats
summaryPPMI_bsl(features(mask_bsl,:),group(mask_bsl),feature_names,PPMI(mask_bsl,:));


%% delete
rmpath('./selected_aux');

end



