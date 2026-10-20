function [mean_results,std_results]=analyze_t_mae_auc(t_mae_auc,feat_names,bstrp)

%% pre-processing
n_features=length(feat_names);
mae_auc=nan(3000,size(t_mae_auc,2),size(t_mae_auc,3));
idx_t=zeros(size(t_mae_auc,3),1);
for i=1:size(t_mae_auc,3)
    sample_mae_auc=t_mae_auc(:,:,i);
    t=ceil(sample_mae_auc(:,1));
    for j=1:size(t_mae_auc,3)
        samples=sum(t==j);
        if(samples>0)
            idx=idx_t(j)+1:idx_t(j)+samples;
            mae_auc(idx,:,j)=sample_mae_auc(t==j,:);
            idx_t(j)=idx_t(j)+samples;
        end
    end
end

%% Results
mean_results=nan(size(mae_auc,3),n_features+3);
std_results=nan(size(mae_auc,3),n_features+3);

for i=1:size(mae_auc,3)
    mae=mae_auc(1:idx_t(i),2:n_features+1,i);
    auc=mae_auc(1:idx_t(i),n_features+2:end,i);
    s_auc=nan(bstrp,1);
    BCA=nan(bstrp,1);
    samples=nan(bstrp,1);
    for n=1:bstrp
        auc_bstrp=[auc(:,1),auc(:,1+n)];
        auc_bstrp=rmmissing(auc_bstrp);
        samples(n)=size(auc_bstrp,1);
        BCA(n)=balanced_accuracy(auc_bstrp(:,1)>=0.5, auc_bstrp(:,2)>=0.5);
        try
            [~,~,~,s_auc(n)]= perfcurve(auc_bstrp(:,1), auc_bstrp(:,2), 1);
        catch
            pause;
        end
    end
    fprintf('Time: %d year, subjects: %.0f, AUC: %.3f, BCA: %.3f, MAE (mean):\t',...
        i,mean(samples),mean(s_auc),mean(BCA));
    for j=1:n_features
        fprintf('%s: %.3f\t',feat_names(j),mean(mae(:,j),'omitmissing'))
    end
    fprintf('\n');

    mean_results(i,:)=[mean(samples),mean(s_auc),mean(BCA),mean(mae,'omitmissing')];
    std_results(i,:)=[std(samples),std(s_auc),std(BCA),std(mae,'omitmissing')];
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BCA = balanced_accuracy(y_true, y_pred)
    % Compute confusion matrix
    C = confusionmat(y_true, y_pred);
    
    % Extract true positive, true negative, false positive, and false negative
    TP = C(1,1);
    TN = C(2,2);
    FP = C(1,2);
    FN = C(2,1);
    
    % Compute sensitivity (true positive rate)
    sensitivity = TP / (TP + FN);
    
    % Compute specificity (true negative rate)
    specificity = TN / (TN + FP);
    
    % Compute balanced accuracy
    BCA = (sensitivity + specificity) / 2;
end