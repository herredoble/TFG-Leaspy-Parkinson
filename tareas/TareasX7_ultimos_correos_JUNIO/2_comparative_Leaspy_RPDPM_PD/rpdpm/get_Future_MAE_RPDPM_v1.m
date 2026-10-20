function t_mae_auc=...
    get_Future_MAE_RPDPM_v1(y,labels,ages,Sigma,A,B,C,D,G,loss_type, fit_type, ...
    optim_options,btstrp,posterior_bayes,future_timepoints)



I=size(ages,1);
K=size(y,1);

t_mae_auc=nan(I,K+btstrp+2,length(future_timepoints));

for j=1:length(future_timepoints)
    mask_future= ~isnan(ages);
    timepoints=nan(I,1);
    ages_aux=ages;
    for i=1:I
        n_samples=sum(mask_future(i,:));
        if((n_samples/2)>=future_timepoints(j))
            mask_future(i,n_samples-future_timepoints(j)+1:n_samples)=0;    
            timepoints(i)=ages(i,n_samples)-ages(i,n_samples-future_timepoints(j));
            ages_aux(i,n_samples-future_timepoints(j)+1:n_samples-1)=nan;
        end
    end
    
    [label_test_fit,ae_test]=get_Future_MAE_RPDPM(y,ages_aux,...
        Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp,mask_future,posterior_bayes);

    mae_subj=nan(I,size(y,1));
    label_subj=nan(I,btstrp+1);
    for i=1:I
        mask_subj_future = mask_future(i,:)==0 & ~isnan(ages_aux(i,:));
        if(sum(mask_subj_future)>0)
            mae_subj(i,:)= mean(squeeze(ae_test(:,i,mask_subj_future,:))','omitmissing');
            label_subj(i,:)=[labels(i,mask_subj_future),squeeze(label_test_fit(i,mask_subj_future,:))'];
        end
    end
    t_mae_auc(:,:,j)=[timepoints,mae_subj,label_subj];
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label_test_fit,ae_test]=...
    get_Future_MAE_RPDPM(y_test,ages_test,Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp,mask_future,posterior_bayes) 
%% Robust Parametric DPM Testing

% Estimating test subject parameters using biomaker values and ages
[K2, I2, J2] = size(y_test); % number of biomarkers, subjects, and visits
y_test_fit = zeros(K2, I2, J2, btstrp); % estimated values of biomarkers
label_test_fit = zeros(I2, J2, btstrp);
Alpha_test = zeros(I2, btstrp);
Beta_test = zeros(I2, btstrp);
for n = 1 : btstrp
    [Alpha_test(:, n), Beta_test(:, n), label_test_fit(:, :, n), y_test_fit(:, :, :, n)] = rpdpm_test_future(y_test, ages_test, ...
        Alpha_test(:, n), Beta_test(:, n), Sigma(:, n), A(:, n), B(:, n), C(:, n), D(:, n), G(:, n), loss_type, fit_type, optim_options,...
        mask_future,posterior_bayes(2,n));

end

% Test set modeling performance
ae_test= nan(K2,I2,J2,btstrp);
for n = 1 : btstrp
    ae_test(:,:,:,n) = abs(y_test - y_test_fit(:, :, :, n)); % absolute errors
end



end
