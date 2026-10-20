function [idx_train,idx_test]= split_data_RPDPM(y,classes,diagnose,ratio_test)
%% Splitting data to training and test subsets

    [~,I,~] = size(y); % number of unique subjects
    num_points = sum(sum(~isnan(y), 3), 1); % number of available points per subject
    idx_test = []; % indices of test subjects
    rng(0); % random number generation seed
    for u = 1 : numel(classes)
        for uu = u : numel(classes)
            idx_u = find(strcmp(diagnose(1, :), classes(u)) & strcmp(diagnose(2, :), classes(uu))); % first and last diagnoses
            med_u = median(num_points(idx_u)); % within class median of available points
            idx_many = find(num_points(idx_u) > med_u); % subjects with many points
            if ratio_test * numel(idx_many) > 1
                idx_rnd = randperm(numel(idx_many), floor(ratio_test * numel(idx_many)));
                idx_test = cat(2, idx_test, idx_u(idx_many(idx_rnd)));
            end
            idx_few = find(num_points(idx_u) <= med_u); % subjects with few points
            if ratio_test * numel(idx_few) > 1
                idx_rnd = randperm(numel(idx_few), floor(ratio_test * numel(idx_few)));
                idx_test = cat(2, idx_test, idx_u(idx_few(idx_rnd)));
            end
        end
    end
    idx_train = find(~ismember(1 : I, idx_test)); % indices of training subjects

end