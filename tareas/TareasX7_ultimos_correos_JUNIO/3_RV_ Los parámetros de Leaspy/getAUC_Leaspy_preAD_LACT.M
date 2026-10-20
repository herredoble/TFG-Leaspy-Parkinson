function getAUC_Leaspy_preAD_LACT
addpath('./aux_grace');
clc;

%% Leaspy parameters
tbl_param_leaspy=readtable('./data/param_leaspy_LACT.xlsx');
RID=str2double(tbl_param_leaspy.ID);
tbl_param_leaspy.RID=RID;


%% ADNI data
n_typeFeat= 2; %CAM
typeFeatures={'ROI','CAM','ROICAM','all'};
type_features=cellstr(typeFeatures{n_typeFeat});
[ADNI,~,~,~,~] = getDataADNIMERGE_v1_2b(type_features);


mask_CN=(ADNI.stMCI==0) | (ADNI.stMCI<0); % CU and SMC
mask_sMCI=ADNI.stMCI==1 & ADNI.Convert==0;
mask=mask_CN | mask_sMCI;
ADNI=ADNI(mask,:);
ADNI=sortrows(ADNI,{'RID','M'});
ADNI(ADNI.DX=='Dementia',:)=[];

%% Get DPS
RID=sort(RID);
J = histcounts(ADNI.RID, RID); % number of available visits per subject
btstrp=1;

DPS = NaN(length(RID), max(J),btstrp);
ages = NaN(length(RID), max(J)); 
labels  = cell(length(RID), max(J));
for i=1:length(RID)
    param_subj=tbl_param_leaspy(tbl_param_leaspy.RID==RID(i),:);
    ADNI_subj= ADNI(ADNI.RID==RID(i),:);
    tij=ADNI_subj.AGE+(ADNI_subj.M)/12;
    DPS(i,1:length(tij),1)= getDPS(tij,param_subj.tau,param_subj.xi)';
    labels(i,1:length(tij))= cellstr(ADNI_subj.DX)';
    ages(i,1:length(tij))=tij; 
end

%% show
% figure(2);clf;
% DPS_1D=DPS(:);
% labels_1D=labels(:);
% mask_nan = ~isnan(DPS_1D);
% DPS_1D=DPS_1D(mask_nan);
% labels_1D =labels_1D(mask_nan);
% histogram(DPS_1D(strcmp(labels_1D,'CN')));hold on;
% histogram(DPS_1D(strcmp(labels_1D,'MCI')));hold off;
% xlabel('DPS');
% legend('CN','MCI');

%% Training Bayesian classifiers per bootstrap

I=length(RID);
ratioTest = 0.2; % proportion of test samples to entire data
idx_test = randperm(I, ceil(ratioTest * I)); % random indices of test samples
idx_train = find(~ismember(1 : I, idx_test)); % random indices of training samples
DPS_train=DPS(idx_train,:);
DPS_test=DPS(idx_test,:,:);
labels_train=labels(idx_train,:);
labels_test=labels(idx_test,:);
ages_train=ages(idx_train,:);
ages_test=ages(idx_test,:);

classes = {'CN','MCI'};
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
[ages_onset,gt_onset,m_test]=get_Onset_MCI_preAD(DPS_test,...
    posterior_bayes,ages_test,labels_test,btstrp);
convert_times=[ages_onset,gt_onset];

%% Train
fprintf('Train \t');
auc_test=get_AUC_RPDPM_v1(DPS_train,posterior_bayes,labels_train,classes,btstrp);
[ages_onset,gt_onset,m_test]=get_Onset_MCI_preAD(DPS_train,...
    posterior_bayes,ages_train,labels_train,btstrp);
convert_times_=[convert_times;ages_onset,gt_onset];
convert_times_=rmmissing(convert_times_);

figure(1);clf;
plot(convert_times_(:,2),convert_times_(:,1),'+');
xlabel('decline cognitive convert age');
ylabel('Leaspy estimation');
c=corrcoef(convert_times_(:,2),convert_times_(:,1));
fprintf('Decline cognitive convert age correlation: %.3f\n',c(1,2));


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://disease-progression-modelling.github.io/pages/notebooks/disease_course_mapping/TP2_leaspy_beginner.html
function DPS = getDPS(tij,tau,xi)
DPS=exp(xi)*(tij-tau);
% DPS=(tij-tau);
end
