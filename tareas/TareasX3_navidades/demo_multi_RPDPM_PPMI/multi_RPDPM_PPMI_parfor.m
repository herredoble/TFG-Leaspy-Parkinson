function multi_RPDPM_PPMI_parfor
clc;
close all;

addpath('./aux_grace');
addpath('./rpdpm'); % add path to the source codes

% delete(gcp('nocreate'));
% parpool(feature('numcores'));


%% Input data file
%addpath('./aux_grace');
n_dataset=1; %PPMIApril20
n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
[PPMI,features,feature_names,~,~] = getDataPPMI_MCI_v1_2(n_typeFeat,n_dataset);
feature_names=string(feature_names)';
feature_names=upper(feature_names);

% PPMI.abeta(PPMI.PATNO==3802)=700;
% PPMI.abeta(PPMI.PATNO==3868)=700;

%feature_tbl_names=string(PPMI.Properties.VariableNames)';
%selectedColums=[2,159,16,166,167,76,103,71,99,112,134];% LNS- SDMT- HVLTRDLY- SFT- MEANCAU- ABETA-
selectedColums=[2,159,16,166,167];
data=PPMI(:,selectedColums);
data.AGE=data.age+years(data.date-data.date_bsl);
% Distinct class labels ordered w.r.t. disease progression
%classes = {'CN', 'MCI', 'AD'};
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
features_set = [1,2,4:10,12,14,15,17,24];
%features_set = [15,9,5,14,6,10];


% Rango de tamaños del vector a investigar
vector_range = 3:5;
%vector_range = 5:6;
%vector_range = 5:7;

head_fich_name ="./data_multi/PPMI_MCI_";
btstrp=10;

%% feature table
feature_names=feature_names(features_set);
% feature_names=erase(feature_names,"RAVLT_");
% feature_names=erase(feature_names,"_");

features_rpdpm = array2table(features(:,features_set),"VariableNames",feature_names);
data_rpdpm=[features_rpdpm,data];

% Distinct class labels ordered w.r.t. disease progression
%classes = {'CN', 'MCI', 'AD'};
classes = {'CN', 'MCI'};

% Ranges of biomarkers
%[data_rpdpm,ranges]=set_ranges_PPMI_LSHSM(data_rpdpm);
[data_rpdpm,ranges]=set_ranges_PPMI(data_rpdpm);

% Proportion of test subjects to entire data
ratio_test = 0.2;





%% Normalized data

[y,ages,labels,diagnose,subjects,idx_train,idx_test]=...
    data_preparation_RPDPM_PPMI(data_rpdpm,ranges,classes,ratio_test);

save('./data/data_PPMI_model_multi','y','ages','labels','diagnose','subjects',...
    'idx_train','idx_test','feature_names');


%% Pool marker sets    

features_set=1:length(features_set);
max_dim=vector_range(end);
features_subset_matrix = nan(10000,max_dim);
next_row=1;
for i = 1 : length(vector_range)
    aux_subsets = nchoosek(features_set,vector_range(i));
    features_subset_aux=nan(size(aux_subsets,1),max_dim);
    features_subset_aux(:,1:vector_range(i))=aux_subsets;
    n_rows=size(features_subset_aux,1);
    features_subset_matrix(next_row:n_rows+next_row-1,:)=features_subset_aux;
    next_row=next_row+n_rows;
end
features_subset_matrix(next_row:end,:)=[];



%%
pc = parcluster('local');
pc.JobStorageLocation = strcat(getenv('SCRATCH'),'/', getenv('SLURM_JOB_ID'));
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

for i = 1 : size(features_subset_matrix,1)
    idx_markers=features_subset_matrix(i,:);
    idx_markers=idx_markers(~isnan(idx_markers));
    %idx_markers=1:6;
    model_feat=strjoin(feature_names(idx_markers), ', ');
    fprintf('Model %d of %s\n', i, model_feat);
    str_feat=strrep(int2str(idx_markers), ' ', '_');
    fich_name=head_fich_name+str_feat;
    if(~exist(fich_name+'.mat','file'))
        try
            [~,~,~,~,~,~,~]=...
                build_rpdpm_PPMI(y(idx_markers,:,:),ages,labels,diagnose,...
                classes,ranges(idx_markers,:),idx_train,idx_test,btstrp,...
                fich_name,feature_names(idx_markers));
        catch
        end
    end
        
    
    fprintf('\n\n');
end
   
    




%% End

rmpath('./aux_grace');
delete(gcp('nocreate'));

end



