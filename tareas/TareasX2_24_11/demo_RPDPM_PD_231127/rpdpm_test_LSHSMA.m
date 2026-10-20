function rpdpm_test_LSHSMA
%% Training, validation, and testing the proposed model for biomarker value prediction and clinical status classification using a simulated data 
% MAE (LNS): 1.538	MAE (SDMTOTAL): 4.994	MAE (HVLTRDLY): 1.548	MAE (SFT): 38.816	MAE (MEANCAU): 0.413	MAE (ABETA): 290.458	
% NMAE = 0.6210 ± 0.0065 
%  AUC = 0.8202 ± 0.0033 
%% https://github.com/Mostafa-Ghazi/CARRNN
restoredefaultpath
close all
clc

addpath('./rpdpm'); % add path to the source codes

%% Input data file
addpath('./aux_grace');

n_dataset=1; %PPMIApril20
n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
[PPMI,~,~,~,~] = getDataPPMI_MCI_v1_2(n_typeFeat,n_dataset);


% PPMI.abeta(PPMI.PATNO==3802)=700;
% PPMI.abeta(PPMI.PATNO==3868)=700;

%feature_tbl_names=string(PPMI.Properties.VariableNames)';
selectedColums=[2,159,16,166,167,76,103,71,99,112,134];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
data=PPMI(:,selectedColums);
data.AGE=data.age+years(data.date-data.date_bsl);
data.DX=categorical(data.cogstate>1,[0,1],{'CN','MCI'});
data.DX=cellstr(data.DX);
data(:,2:5)=[];
data.Properties.VariableNames{'PATNO'} = 'RID';
data.Properties.VariableNames{'lns'} = 'LNS';
data.Properties.VariableNames{'sft'} = 'SFT';
data.Properties.VariableNames{'abeta'} = 'ABETA';
data.Properties.VariableNames{'mean_caudate'} = 'MEANCAU';
data=data(:,[1,7,4,2,6,3,5,8,9]);

% Ranges of biomarkers
data.SDMTOTAL(data.SDMTOTAL>80)=80;
data.ABETA(data.ABETA>2000)=2000;
data.SFT(data.SFT>100)=100;
%ranges = [0 20; 0 80;0 12;7 13;0.1650 3.7050;200 2000];
ranges = [200 2000;0 12;0 20;0.1650 3.7050; 0 80;7 100];

colum_names = string(data.Properties.VariableNames);
feat_names=colum_names(2:7);

%% Incresead features
selectedFeatures=2:7;
features=table2array(data(:,selectedFeatures));
change_sign_feat=true(length(selectedFeatures),1);  
for i=1:length(change_sign_feat)
    if(change_sign_feat(i)==true)
        features(:,i)=-features(:,i);
    end
end
%% Transforming
group=strcmp(data.DX,'MCI');
[~,percent,point_feat]=percentil_interp(features,group);
%[f2,percent,point_feat]=percentil_interp(features,[],percent,point_feat);

%% Parameters

% Distinct class labels ordered w.r.t. disease progression
classes = {'CN', 'MCI'};

% Optimization parameters
loss_type = 'logistic'; % robust estimation loss function type
fit_type = 'proposed'; % fitting function type
btstrp = 10; % number of bootstraps

%% Data Preparation
% Extracting desired fields from the longitudinal data
subjects_vec = data{:, 'RID'}; % subject IDs
labels_vec = data{:, 'DX'}; % visiting status of subjects
ages_vec = data{:, 'AGE'}; % visiting ages of subjects
biomarkers = setdiff(data.Properties.VariableNames, {'RID', 'DX', 'AGE'}); % biomarker names
y_vec = data{:, ismember(data.Properties.VariableNames, biomarkers)}; % biomarker values


% Arranging data in a desired form of 2D/3D arrays
subjects = unique(subjects_vec, 'sorted'); % unique subjects
I = numel(subjects); % number of unique subjects
K = numel(biomarkers); % number of biomarkers
J = histcounts(subjects_vec, I); % number of available visits per subject
y = NaN(K, I, max(J)); % biomarker values
ages = NaN(I, max(J)); % visiting ages of subjects
labels = cell(I, max(J)); % visiting status of subjects
for i = 1 : I
    idx_i = find(subjects_vec == subjects(i));
    [val_sort, idx_sort] = sort(ages_vec(idx_i), 'ascend');
    ages(i, 1 : numel(idx_i)) = val_sort;
    y(:, i, 1 : numel(idx_i)) = y_vec(idx_i(idx_sort), :)';
    labels(i, 1 : numel(idx_i)) = labels_vec(idx_i(idx_sort));
end

% Finding clinical status of subjects in their first and last available visits
% diagnose =  clinical_group_first_end(labels);


%% Robust Parametric DPM Training (model)
% Optimization options
optim_options = optimoptions('lsqnonlin', 'Display', 'none', 'Jacobian', 'off', 'DerivativeCheck', 'off', 'MaxIter', 1e+3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Algorithm', 'trust-region-reflective');

load('./data/workspace_RPDPM_PD_LSHSMA_10_bstrp.mat', 'Sigma', 'A', 'B', 'C', 'D', 'G','posterior_bayes' );
warning('off');

% [Sigma, A, B, C, D, G, DPS_train, Idx_train, idx_valid]=...
%     build_model_RPDPM(y_train,ages_train,labels_train,diagnose_train,...
%     classes,ranges,loss_type, fit_type, optim_options,btstrp,L_max,L_min);
% 
% 
%% Robust Parametric DPM Testing
n_subjs=length(subjects);
%n_subjs=5;

y_test=y(:,1:n_subjs,:);
ages_test=ages(1:n_subjs,:);
RID_test=subjects(1:n_subjs);

years_before=1;
years_after=5;
resolution_year=0.2;

timeline=get_timeline_tray(ages_test,years_before,years_after,resolution_year);


[DPS_test,DPS_tray,y_test_fit,y_tray_fit]=get_test_tray_RPDPM(y_test,feat_names,ages_test,timeline,...
    Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp);

% tSample=buildTableSample(RID_test,y_test,y_test_fit,ages_test,DPS_test,string(biomarkers));
% tTrayec=buildTableSample(RID_test,[],y_tray_fit,timeline,DPS_tray,string(biomarkers));
tSample=buildTableSample(RID_test,y_test,y_test_fit,ages_test,DPS_test,feat_names);
tTrayec=buildTableSample(RID_test,[],y_tray_fit,timeline,DPS_tray,feat_names);

%% Robust Parametric DPM Classification Testing
labels_test=labels(1:n_subjs,:)';
[predict_scores_sample,labels_test_1D]=prob_RPDPM(tSample,posterior_bayes,labels_test,classes,btstrp);
[predict_scores_trayec,~]=prob_RPDPM(tTrayec,posterior_bayes,[],classes,btstrp);


tSample.DX=labels_test_1D;
for n=1:btstrp
    probMCI=predict_scores_sample(2,:,n)';
    tSample=[tSample,table(probMCI,'VariableNames',string(sprintf('probMCI_%d',n)))];    
    probMCI=predict_scores_trayec(2,:,n)';
    tTrayec=[tTrayec,table(probMCI,'VariableNames',string(sprintf('probMCI_%d',n)))];    

end

% [predict_scores_trayec,~]=prob_RPDPM(tTrayectories.DPS,posterior_bayes,[],classes,btstrp);
% tTrayectories.probCU=predict_scores_trayec(1,:)';
% tTrayectories.probMCI=predict_scores_trayec(2,:)';

%% Cognitive reserve
tau=nan(length(RID_test),btstrp);
for i=1:length(RID_test)
    mask_subj=tTrayec.RID==RID_test(i);
    idx_subj_first=find(mask_subj==1,1);
    %idx_subj_last=find(mask_subj==1,1,'last');
    for n =1:btstrp
        probMCI=table2array(tTrayec(:, string(sprintf('probMCI_%d',n))));
        idx=find(probMCI(mask_subj)>0.5,1);
        if(isempty(idx))%CU even tmax+5
            tau(i,n)=inf;%tTrayec.Age(idx_subj_last);
        elseif(idx==1)%MCI even 5 years before
             tau(i,n)=tTrayec.Age(idx_subj_first);
        else
            tau(i,n)=(tTrayec.Age(idx+idx_subj_first-1)+tTrayec.Age(idx+idx_subj_first-2))/2;
        end
    end
end

end




