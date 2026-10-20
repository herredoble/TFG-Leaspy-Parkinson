function posterior_bayes = get_AUC_Bayes (DPS_train,labels_train,classes,btstrp)

% Training Bayesian classifiers per bootstrap
UC = numel(classes); % number of unique class labels
dist_models = repmat({'kernel'}, 1, UC); % distribution models for likelihoods
posterior_bayes = cell(UC, btstrp);
likelihood_bayes = cell(UC, btstrp);
evidence_bayes = cell(1, btstrp);
prior_bayes = zeros(UC, btstrp);
for n = 1 : btstrp
    [posterior_bayes(:, n), likelihood_bayes(:, n), prior_bayes(:, n), evidence_bayes{n}]...
        = bayes_train(DPS_train(:, :, n), labels_train, classes, dist_models);
end

end
