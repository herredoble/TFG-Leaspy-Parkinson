function build_rpdpm_LSHSMA
%% Training, validation, and testing the proposed model for biomarker value prediction and clinical status classification using a simulated data 
% Training BIC = 3718.5672 ± 92.6179 
% Validation NMAE = 0.6318 ± 0.0162 
% Test 	MAE (LNS): 1.453	MAE (SDMTOTAL): 4.827	MAE (HVLTRDLY): 1.560	MAE (SFT): NaN	MAE (MEANCAU): 0.379	MAE (ABETA): 284.774	
% NMAE = 0.5972 ± 0.0069 
% Validation 	MAE (LNS): 1.524	MAE (SDMTOTAL): 4.718	MAE (HVLTRDLY): 1.622	MAE (SFT): 2.034	MAE (MEANCAU): 0.401	MAE (ABETA): 303.865	
% NMAE = 0.6208 ± 0.0063 
% Train 	MAE (LNS): 1.564	MAE (SDMTOTAL): 5.059	MAE (HVLTRDLY): 1.545	MAE (SFT): 2.488	MAE (MEANCAU): 0.417	MAE (ABETA): 293.367	
% NMAE = 0.6271 ± 0.0064 
% Test 	 AUC = 0.8435 ± 0.0087 
% Validation 	 AUC = 0.8740 ± 0.0030 
% Train 	 AUC = 0.8183 ± 0.0032 

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
data.SDMTOTAL(data.SDMTOTAL>80)=80;
data.ABETA(data.ABETA>2000)=2000;
data.SFT(data.SFT>100)=100;

colum_names = string(data.Properties.VariableNames);
feat_names=colum_names(2:7);



%% Parameters
% Ranges of biomarkers
ranges = [200 2000;0 12;0 20;0.1650 3.7050; 0 80;7 100];

% Proportion of test subjects to entire data
ratio_test = 0.2;

% Distinct class labels ordered w.r.t. disease progression
%classes = {'CN', 'MCI', 'AD'};
classes = {'CN', 'MCI'};

% Optimization parameters
loss_type = 'logistic'; % robust estimation loss function type
fit_type = 'proposed'; % fitting function type
L_max = 50; % maximum number of alternating iterations
L_min = 10; % minimum number of alternating iterations
btstrp = 10; % number of bootstraps

%% Data Preparation

% Extracting desired fields from the longitudinal data
subjects_vec = data{:, 'RID'}; % subject IDs
labels_vec = data{:, 'DX'}; % visiting status of subjects
ages_vec = data{:, 'AGE'}; % visiting ages of subjects
biomarkers = setdiff(data.Properties.VariableNames, {'RID', 'DX', 'AGE'}); % biomarker names
y_vec = data{:, ismember(data.Properties.VariableNames, biomarkers)}; % biomarker values

% Removing data lacking visiting age information and vice versa
y_vec(isnan(ages_vec), :) = NaN;
ages_vec(~sum(~isnan(y_vec), 2)) = [];
labels_vec(~sum(~isnan(y_vec), 2)) = [];
subjects_vec(~sum(~isnan(y_vec), 2)) = [];
y_vec(~sum(~isnan(y_vec), 2), :) = [];

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
diagnose =  clinical_group_first_end(labels);

% Filtering data by rejecting outliers
for k = 1 : K
    idx_ij = find(~isnan(y(k, :)));
    y_k = y(k, idx_ij);
    idx_out = find(y_k < ranges(k, 1) | y_k > ranges(k, 2));
    y(k, idx_ij(idx_out)) = NaN; % reject outliers
end

% Removing subjects with a few number of distinct visits
idx_rmv = find(sum(~isnan(nanmean(y, 1)), 3) < 2 | sum(~cellfun(@isempty, labels), 2)' < 2);
y(:, idx_rmv, :) = [];
ages(idx_rmv, :) = [];
labels(idx_rmv, :) = [];
diagnose(:, idx_rmv) = [];
subjects(idx_rmv, :) = [];
%I = numel(subjects); % number of unique subjects


% Removing subjects with MCI first and CN last
idx_rmv = find(strcmp(diagnose(1, :), classes(2)) & strcmp(diagnose(2, :), classes(1)));
y(:, idx_rmv, :) = [];
ages(idx_rmv, :) = [];
labels(idx_rmv, :) = [];
diagnose(:, idx_rmv) = [];
subjects(idx_rmv, :) = [];


%% Splitting data to training and test subsets
[idx_train,idx_test]= split_data_RPDPM(y,classes,diagnose,ratio_test);

y_train = y(:, idx_train, :);
ages_train = ages(idx_train, :);
labels_train = labels(idx_train, :);
diagnose_train = diagnose(:, idx_train);
subjects_train = subjects(idx_train, :);

y_test = y(:, idx_test, :);
ages_test = ages(idx_test, :);
labels_test = labels(idx_test, :);
diagnose_test = diagnose(:, idx_test);
subjects_test = subjects(idx_test, :);


%% Robust Parametric DPM Training (model)
% Optimization options
optim_options = optimoptions('lsqnonlin', 'Display', 'none', 'Jacobian', 'off', 'DerivativeCheck', 'off', 'MaxIter', 1e+3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Algorithm', 'trust-region-reflective');


[Sigma, A, B, C, D, G, DPS_train, Idx_train, idx_valid]=...
    build_model_RPDPM(y_train,ages_train,labels_train,diagnose_train,...
    classes,ranges,loss_type, fit_type, optim_options,btstrp,L_max,L_min);



% Visualizing the DPM results
% plot_marker_order(C, biomarkers);
% plot_dps_pdf(DPS_train, C, labels_train, classes, biomarkers, Idx_train);
% plot_trajectory(A, B, C, D, G, DPS_train, y_train, labels_train, classes, biomarkers, Idx_train, fit_type);

% plot_marker_order(C, feat_names);
% plot_trajectory(A, B, C, D, G, DPS_train, y_train, labels_train, classes, feat_names, Idx_train, fit_type);

%% Robust Parametric DPM Testing
fprintf('Test \t');
DPS_test=get_MAE_RPDPM(y_test,feat_names,ages_test,Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp);

%% Robust Parametric DPM Validation
y_val=y_train(:, idx_valid, :);
ages_val = ages_train(idx_valid, :);
fprintf('Validation \t');
DPS_val=get_MAE_RPDPM(y_val,feat_names,ages_val,Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp);

%% Robust Parametric DPM train
fprintf('Train \t');
DPS_train_2=get_MAE_RPDPM(y_train,feat_names,ages_train,Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp);

%% Robust Parametric DPM Classification Training
% Training Bayesian classifiers per bootstrap
UC = numel(classes); % number of unique class labels
dist_models = repmat({'kernel'}, 1, UC); % distribution models for likelihoods
posterior_bayes = cell(UC, btstrp);
likelihood_bayes = cell(UC, btstrp);
evidence_bayes = cell(1, btstrp);
prior_bayes = zeros(UC, btstrp);
for n = 1 : btstrp
    [posterior_bayes(:, n), likelihood_bayes(:, n), prior_bayes(:, n), evidence_bayes{n}] = bayes_train(DPS_train(:, :, n), labels_train(Idx_train(:, n), :), classes, dist_models);
end

%% Robust Parametric DPM Classification Testing
fprintf('Test \t');
get_AUC_RPDPM(DPS_test,posterior_bayes,labels_test,classes,btstrp);

%% Robust Parametric DPM Classification Validation
labels_val = labels_train(idx_valid, :);
fprintf('Validation \t');
get_AUC_RPDPM(DPS_val,posterior_bayes,labels_val,classes,btstrp);

%% Robust Parametric DPM Classification Validation
fprintf('Train \t');
get_AUC_RPDPM(DPS_train_2,posterior_bayes,labels_train,classes,btstrp);

end




