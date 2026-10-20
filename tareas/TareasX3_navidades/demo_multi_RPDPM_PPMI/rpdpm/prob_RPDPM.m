function [predict_scores_test,labels_test_1D]=prob_RPDPM(tSample,posterior_bayes,labels_test,classes,btstrp)
%% Robust Parametric DPM Classification Testing
if(~isempty(labels_test))
    labels_test_1D=cell(size(tSample,1),1);
    I=length(unique(tSample.RID));
    r=1;
    for i=1:I
        %n_visits_y=sum(tSample.RID==RID_test(i));
        vacias = cellfun(@isempty, labels_test(:,i));
        n_visits_l=sum(~vacias);
        %if(n_visits_y~=n_visits_l)
        %    fprintf('Visits are diferents in %d\n',RID(i));
        %end
        labels_test_1D(r:r+n_visits_l-1)=labels_test(vacias==0,i);
        r=r+n_visits_l;
    end
else
    labels_test_1D=[];   
end

% Classification performance analysis for test samples
UC = numel(classes); % number of unique class label
n_visits=size(tSample,1);
predict_scores_test = nan(UC, n_visits, btstrp); % posterior probabilities for test scores
auc_test = zeros(1, btstrp); % test multiclass area under the curve per bootstrap
for n = 1 : btstrp
    
    dps_test_n =  table2array(tSample(:, string(sprintf('dps_%d',n))));
    for cc = 1 : UC
        predict_scores_test(cc, :, n) = posterior_bayes{cc, n}(dps_test_n);
    end
    if(~isempty(labels_test))
        auc_test(n) = multiclass_auc(predict_scores_test(:, :, n)', ...
            labels_test_1D, classes);       
    end
end
if(~isempty(labels_test))
    auc_test_avg = mean(auc_test); % average multiclass area under the curve
    auc_test_std = std(auc_test); % standard deviation of multiclass areas under the curves
    fprintf(' AUC = %4.4f \x00B1 %4.4f \n', [auc_test_avg, auc_test_std]); % display test classification performance
end

end
