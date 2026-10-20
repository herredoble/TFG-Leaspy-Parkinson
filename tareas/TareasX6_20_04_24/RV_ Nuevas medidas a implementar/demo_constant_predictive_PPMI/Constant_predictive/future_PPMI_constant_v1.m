function future_PPMI_constant_v1
clc
addpath('./aux_grace');
addpath('./rpdpm'); % add path to the source codes

%% Input data file
n_dataset=1; %PPMIApril20
n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
[PPMI,features,feature_names,~,~] = getDataPPMI_MCI_v1_2(n_typeFeat,n_dataset);
feature_names=string(feature_names)';
feature_names=upper(feature_names);


%feature_tbl_names=string(PPMI.Properties.VariableNames)';
%selectedColums=[2,159,16,166,167,76,103,71,99,112,134];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
selectedColums=[2,159,16,166,167];
data=PPMI(:,selectedColums);
data.AGE=data.age+years(data.date-data.date_bsl);
classes = {'CN', 'MCI'};

data.DX=categorical(data.cogstate>1,[0,1],classes);
data.DX=cellstr(data.DX);
data(:,2:5)=[];
data.Properties.VariableNames{'PATNO'} = 'RID';

% Markers:
%     1. MOCA(-)
%     2. GDS (+)
%     3. UPSIT (-)
%     4. REM (+)
%     5. LNS (-)
%     6. SDMT (-)
%     7. JLO (-)
%     8. INMEDIATE (-)
%     9. HVLTRDLY (-)
%     10. SFT (-)
%     11. ADL (-)
%     12. SCOPA (+)
%     13. MEANPUT (-)
%     14. MEANCAU (-)
%     15. ABETA (-)
%     16. ASYN (-)
%     17. TAU (+)
%     18. PTAU (+)
%     19. pTauAB (+)
%     20. UPDRSTOT
%     21. UPDRSTOTON
%     22. UPDRS1
%     23. UPDRS2
%     24. UPDRS3 (+)
%     25. UPDRS3ON

%features_set = [1,2,4:10,14,19,24];
%features_set = [1,2,4:10,14,15];
features_set = [1,2,4:10,12,14,15,17,24]; %Dic23
%features_set = [15,9,5,14,6,10]; %LSHSM

%% feature table
feature_names=feature_names(features_set);
features_rpdpm = array2table(features(:,features_set),"VariableNames",feature_names);
data_rpdpm=[features_rpdpm,data];


% Ranges of biomarkers
[data_rpdpm,ranges]=set_ranges_PPMI(data_rpdpm);%Dic23

% Proportion of test subjects to entire data
ratio_test = 0.2;

[y,ages,labels,~,~,idx_train,idx_test]=...
    data_preparation_RPDPM_PPMI(data_rpdpm,ranges,classes,ratio_test);

%% Split data
labels_num=nan(size(labels));
labels_num(strcmp(labels,'CN'))=0;
labels_num(strcmp(labels,'MCI'))=1;

labels_train=labels_num(idx_train,:);
labels_test=labels_num(idx_test,:);

y_train=y(:,idx_train,:);
y_test=y(:,idx_test,:);

ages_train=ages(idx_train,:);
ages_test=ages(idx_test,:);

%subjects_train=subjects(idx_train);
%subjects_test=subjects(idx_test);

%% Model
selectedFeatures=[6;7;9];%JIS
y_test=y_test(selectedFeatures,:,:);
y_train=y_train(selectedFeatures,:,:);
y=y(selectedFeatures,:,:);
feature_names=feature_names(selectedFeatures);
%% Constant prediction
bstrp=1;
future_timepoints=1:3;
fprintf('Test \n');
t_mae_auc_test=get_Future_MAE_Constant_v1(y_test,labels_test,ages_test,future_timepoints);
[~,~]=analyze_t_mae_auc(t_mae_auc_test,feature_names,bstrp);
fprintf('Train \n');
t_mae_auc_train=get_Future_MAE_Constant_v1(y_train,labels_train,ages_train,future_timepoints);
[~,~]=analyze_t_mae_auc(t_mae_auc_train,feature_names,bstrp);
fprintf('Total \n');
t_mae_auc_total=get_Future_MAE_Constant_v1(y,labels_num,ages,future_timepoints);
[~,~]=analyze_t_mae_auc(t_mae_auc_total,feature_names,bstrp);

end



