function diagnose =  clinical_group_first_end(labels)

    I=size(labels,1);
    
    % Finding clinical status of subjects in their first and last available visits
    idx_labels = ~cellfun(@isempty, labels); % logical indices of available labels
    labels_exist = labels(sum(idx_labels, 2) > 0, :); % labels having at least one label per subject
    idx_labels_exist = idx_labels(sum(idx_labels, 2) > 0, :); % logical indices of available labels having at least one label per subject
    idx_labels_first = arrayfun(@(x) find(idx_labels_exist(x, :), 1, 'first'), 1 : size(idx_labels_exist, 1)); % index of first label per subject
    idx_labels_first = sub2ind(size(idx_labels_exist), 1 : size(idx_labels_exist, 1), idx_labels_first); % linear array indices of first available labels
    idx_labels_last = arrayfun(@(x) find(idx_labels_exist(x, :), 1, 'last'), 1 : size(idx_labels_exist, 1)); % index of last label per subject
    idx_labels_last = sub2ind(size(idx_labels_exist), 1 : size(idx_labels_exist, 1), idx_labels_last); % linear array indices of last available labels
    diagnose = cell(2, I); % first and last diagnosis of subjects
    diagnose(:, sum(idx_labels, 2) > 0) = [labels_exist(idx_labels_first); labels_exist(idx_labels_last)];

end