function summaryPPMI_PD_features(PDtrue_HCfalse,n_Population)
%function summaryPPMI_PD_features(n_Population,n_cog_criteria)
% function summaryPPMI_PD_features(n_Population,n_cog_criteria, bAmyloid)
PDtrue_HCfalse=1;
n_Population=1;

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
    [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_MCI_v1(n_typeFeat,...
        n_cohort,PDtrue_HCfalse);
elseif(strcmp(typeDataset,'PPMISep20_PD_HC'))
    [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_HC_PD(n_typeFeat,...
        n_cohort);
else
end
    
% Cognitive criteria
if(strcmp(varCogCriteria,'SITE'))
    %PPMI.MCI=double(PPMI.cogstate_end>1);
    %group=categorical(PPMI.MCI,[1 0],{'PD-MCI+' 'PD-MCI-'});
    PPMI.MCI=PPMI.cogstate_end;
    PPMI.convertTime=PPMI.convertTimeCogstate;
    if(strcmp(typeDataset,'PPMISep20_PD'))
        if(PDtrue_HCfalse)
            group=categorical(PPMI.MCI,[1 2 3],{'PD-CU' 'PD-MCI' 'PD-Dementia'});
        else
            group=categorical(PPMI.MCI,[1 2 3],{'HC-CU' 'HC-MCI' 'HC-Dementia'});
        end
    elseif(strcmp(typeDataset,'PPMISep20_PD_HC'))
        group=categorical(PPMI.APPRDX,[1 2],{'PD' 'HC'});
    else
    end
% elseif(strcmp(varCogCriteria,'MCI_test'))
%     PPMI.MCI=double(PPMI.MCI_testscores_end);
%     PPMI.convertTime=PPMI.convertTimeMCI_test;
%     group=categorical(PPMI.MCI,[1 0],{'PD-MCI+' 'PD-MCI-'});
% elseif(strcmp(varCogCriteria,'UPDRS1'))
%     PPMI.MCI=double(PPMI.np1cog_end>1);
%     PPMI.convertTime=PPMI.convertTimeNp1cog;
%     group=categorical(PPMI.MCI,[1 0],{'PD-MCI+' 'PD-MCI-'});
else
end

    
        
% elseif(strcmp(typeDataset,'PPMISep20_PD_HC'))
%     n_typeFeat=5; %typeFeatures={'NM','DAT','CSF','MOTOR','all'};
%     n_cohort=1;
%     [PPMI,features,feature_names,~,features_SNan] = getDataPPMI_HC_PD(n_typeFeat,n_cohort);
%     if(bAmyloid)
%         PPMI.Apos = getApos_PPMI(PPMI);
%         mask_Amyloid=~isnan(PPMI.Apos);
%         PPMI=PPMI(mask_Amyloid,:);
%         features=features(mask_Amyloid,:);
%         features_SNan=features_SNan(mask_Amyloid,:);
%         grp=zeros(size(PPMI,1),1);
%         grp(PPMI.APPRDX ==1 & PPMI.Apos==0)=1;
%         grp(PPMI.APPRDX ==1 & PPMI.Apos==1)=2;
%         grp(PPMI.APPRDX ==2 & PPMI.Apos==0)=3;
%         grp(PPMI.APPRDX ==2 & PPMI.Apos==1)=4;
%         group=categorical(grp,[1 2 3 4],{'PD-A-','PD-A+','HC-A-','HC-A+'});
%     else
%         group=categorical(PPMI.APPRDX,[1 2],{'PD' 'HC'});
%     end

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

%% boxplot & LOO
% selectedFeatures=[1,3,4,5,7,8,9,11,13,18];
% numFeatures=length(selectedFeatures);
% feature_bsl=features(mask_bsl,selectedFeatures);
% group_bsl=double(PPMI.MCI(mask_bsl)>1);
% feature_names=feature_names(selectedFeatures);
% %figure(1);clf;
% for i=1:numFeatures
%     %subplot(1,numFeatures,i);
%     mask_outcome=~isnan(feature_bsl(:,i));
%     %boxplot_withPoints(feature_bsl(mask_outcome,i),group_bsl(mask_outcome),...
%     %    feature_names(i));
%     %fprintf('Scores %s: %.2f %.2f %.2f %.3f\n',...
%     %    feature_names(i),LOO(feature_bsl(mask_outcome,i),group_bsl(mask_outcome)));
%     fprintf('Scores %s: %.2f %.2f %.2f %.3f; cutoff: %.2f (%.2f)\n',...
%         feature_names(i),LOO_v1_1(feature_bsl(mask_outcome,i),group_bsl(mask_outcome)));
%     
% end



%% delete
rmpath('./selected_aux');

end


function boxplot_withPoints(measure,group,name)
boxplot(measure,group,'Notch','on','Labels',{'PD-MCI-','PD-MCI+'});
title(name)
hold on;
group=group+1;
f1=scatter(group+((rand(length(group),1)-.5)/5),measure,'k','filled');
f1.MarkerFaceAlpha = 0.05;
hold off;
end

function quality = LOO(X,diagnose)
typeClass ='linear';
% t = templateDiscriminant();
% t = templateNaiveBayes();
% t = templateSVM('Standardize',1,'KernelFunction','gaussian');
N=size(X,1);
labels=zeros(N,1);
scores=zeros(N,2);

for i=1:N
    X_LOO=X;
    X_LOO(i,:)=[];
    Y_LOO=diagnose;
    Y_LOO(i)=[];
%     model_tr = fitcecoc(X_LOO,Y_LOO,'Learners',t);
%     labels(i)=predict(model_tr,X(i,:));
%     [labels(i),~,~,scores(i,:)]=predict(model_tr,X(i,:));
    [labels(i),~,scores(i,:)] = classify(X(i,:),X_LOO,Y_LOO,typeClass);
end
[C_all,~] = confusionmat(diagnose,labels);

TN=C_all(1,1);
FN=C_all(2,1);
FP=C_all(1,2);
TP=C_all(2,2);

sensitivy=TP/(TP+FN);
specificity=TN/(TN+FP);
PPV = TP/(TP+FP);
NPV= TN/(TN+FN);
ACC= (TN+TP)/(TN+TP+FN+FP);

% LRplus=sensitivy/(1-specificity);
% LRminor=(1-sensitivy)/specificity;

class_labels=unique(diagnose);
[~,~,~,auc]=perfcurve(diagnose,scores(:,2),class_labels(2));

% quality=[round([sensitivy;specificity;PPV;NPV;ACC]*100);LRplus;LRminor;auc];
% quality=[round([sensitivy;specificity;PPV;NPV;ACC]*100);auc];
quality=[([sensitivy;specificity;ACC]*100);auc];

end

function quality = LOO_v1_1(X,diagnose)
typeClass ='linear';
% t = templateDiscriminant();
% t = templateNaiveBayes();
% t = templateSVM('Standardize',1,'KernelFunction','gaussian');
N=size(X,1);
labels=zeros(N,1);
scores=zeros(N,2);
cutoff=zeros(N,1);

for i=1:N
    X_LOO=X;
    X_LOO(i,:)=[];
    Y_LOO=diagnose;
    Y_LOO(i)=[];
%     model_tr = fitcecoc(X_LOO,Y_LOO,'Learners',t);
%     labels(i)=predict(model_tr,X(i,:));
%     [labels(i),~,~,scores(i,:)]=predict(model_tr,X(i,:));
    [labels(i),~,scores(i,:),~,coeff] = classify(X(i,:),X_LOO,Y_LOO,typeClass);
    cutoff(i)=coeff(1,2).const/(-coeff(1,2).linear);
end
[C_all,~] = confusionmat(diagnose,labels);

TN=C_all(1,1);
FN=C_all(2,1);
FP=C_all(1,2);
TP=C_all(2,2);

sensitivy=TP/(TP+FN);
specificity=TN/(TN+FP);
PPV = TP/(TP+FP);
NPV= TN/(TN+FN);
ACC= (TN+TP)/(TN+TP+FN+FP);

% LRplus=sensitivy/(1-specificity);
% LRminor=(1-sensitivy)/specificity;

class_labels=unique(diagnose);
[~,~,~,auc]=perfcurve(diagnose,scores(:,2),class_labels(2));

% quality=[round([sensitivy;specificity;PPV;NPV;ACC]*100);LRplus;LRminor;auc];
% quality=[round([sensitivy;specificity;PPV;NPV;ACC]*100);auc];
quality=[([sensitivy;specificity;ACC]*100);auc;mean(cutoff);std(cutoff)];

end



