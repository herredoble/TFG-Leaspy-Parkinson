function [y,ages,labels,diagnose,subjects,idx_train,idx_test]=data_preparation_RPDPM_PPMI(data,ranges,classes,ratio_test)
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




end