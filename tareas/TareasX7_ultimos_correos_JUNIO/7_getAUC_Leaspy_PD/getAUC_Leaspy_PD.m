function getAUC_Leaspy_PD
addpath('./aux_grace');
addpath('./rpdpm');
clc;

%% RPDPM
% MOCA (1), GDS(2), REM(4) y LNS(5)
%"1__2__4__5"	CU: (97.66	98.04)	MCI: (25	40.98)	Age:(0.965	0.975)
%Reserve:(0.770	0.2320) AUC:(0.841	0.833)

%% Leaspy
% Test 	 AUC = 0.7477 ± 0.0000 
% Percentage PD-CU: 100.0 (0.00)
% Percentage PD-MCI/PDD: 0.0 (0.00)
% Correlation age conversion: NaN (NaN)
% Correlation cognition reserve: NaN (NaN)
% Train 	 AUC = 0.7625 ± 0.0000 
% Percentage PD-CU: 100.0 (0.00)
% Percentage PD-MCI/PDD: 0.0 (0.00)
% Correlation age conversion: NaN (NaN)
% Correlation cognition reserve: NaN (NaN)
% Dementia convert age correlation: NaN

%% Leaspy parameters
tTau_Leaspy_test=readtable('./data/IPs_media_test (tau_xi).xlsx');
tTau_Leaspy_train=readtable('./data/IPs_media_train (tau_xi).xlsx');
tbl_param_leaspy=[tTau_Leaspy_test;tTau_Leaspy_train];
tbl_param_leaspy.RID=str2double(tbl_param_leaspy.ID);

%% ADNI data

n_dataset=1; %PPMIApril20
n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
[PPMI,~,~,~,~] = getDataPPMI_MCI_v1_2(n_typeFeat,n_dataset);

%% Get DPS
RID=sort(tbl_param_leaspy.RID);
J = histcounts(PPMI.PATNO, RID); % number of available visits per subject
btstrp=1;

DPS = NaN(length(RID), max(J),btstrp);
ages = NaN(length(RID), max(J)); 
labels  = cell(length(RID), max(J));
for i=1:length(RID)
    param_subj=tbl_param_leaspy(tbl_param_leaspy.RID==RID(i),:);
    PPMI_subj= PPMI(PPMI.PATNO==RID(i),:);
    tij=PPMI_subj.age+PPMI_subj.YEAR;
    DPS(i,1:length(tij),1)= getDPS(tij,param_subj.tau,param_subj.xi)';
    DX=categorical(PPMI_subj.cogstate>1,[0 1],{'CU','MCI'}); 
    labels(i,1:length(tij))= cellstr(DX)';
    ages(i,1:length(tij))=tij; 
end

%% show
figure(1);clf;
DPS_1D=DPS(:);
labels_1D=labels(:);
mask_nan = ~isnan(DPS_1D);
DPS_1D=DPS_1D(mask_nan);
labels_1D =labels_1D(mask_nan);
histogram(DPS_1D(strcmp(labels_1D,'CU')));hold on;
histogram(DPS_1D(strcmp(labels_1D,'MCI')));hold off;
xlabel('DPS');
legend('CU','MCI');

%% Training Bayesian classifiers per bootstrap

idx_test=1:size(tTau_Leaspy_test,1);
idx_train=idx_test(end)+1:size(tbl_param_leaspy,1);
DPS_train=DPS(idx_train,:);
DPS_test=DPS(idx_test,:,:);
labels_train=labels(idx_train,:);
labels_test=labels(idx_test,:);
ages_train=ages(idx_train,:);
ages_test=ages(idx_test,:);

classes = {'CU','MCI'};
UC = numel(classes); % number of unique class labels
dist_models = repmat({'kernel'}, 1, UC); % distribution models for likelihoods
posterior_bayes = cell(UC, btstrp);
likelihood_bayes = cell(UC, btstrp);
evidence_bayes = cell(1, btstrp);
prior_bayes = zeros(UC, btstrp);
for n = 1 : btstrp
    [posterior_bayes(:, n), likelihood_bayes(:, n), prior_bayes(:, n), evidence_bayes{n}] =...
        bayes_train(DPS_train(:, :, n), labels_train(:, :), classes, dist_models);
end
%% Test
fprintf('Test \t');
auc_test=get_AUC_RPDPM_v1(DPS_test,posterior_bayes,labels_test,classes,btstrp);
[ages_onset,gt_onset,m_test]=get_Onset_MCI_PPMI(DPS_test,...
    posterior_bayes,ages_test,labels_test,btstrp);
convert_times=[ages_onset,gt_onset];

%% Train
fprintf('Train \t');
auc_test=get_AUC_RPDPM_v1(DPS_train,posterior_bayes,labels_train,classes,btstrp);
[ages_onset,gt_onset,m_test]=get_Onset_MCI_PPMI(DPS_train,...
    posterior_bayes,ages_train,labels_train,btstrp);
convert_times_=[convert_times;ages_onset,gt_onset];
convert_times_=rmmissing(convert_times_);

figure(2);clf;
plot(convert_times_(:,2),convert_times_(:,1),'+');
xlabel('Dementia convert age');
ylabel('Leaspy estimation');
c=corrcoef(convert_times_(:,2),convert_times_(:,1));
fprintf('Dementia convert age correlation: %.3f\n',c(1,2));


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://disease-progression-modelling.github.io/pages/notebooks/disease_course_mapping/TP2_leaspy_beginner.html
function DPS = getDPS(tij,tau,xi)
DPS=exp(xi)*(tij-tau);
% DPS=(tij-tau);
end
