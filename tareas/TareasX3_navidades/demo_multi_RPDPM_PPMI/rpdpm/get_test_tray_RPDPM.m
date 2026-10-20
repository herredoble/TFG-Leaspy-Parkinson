function [DPS_test,DPS_tray,y_test_fit,y_tray_fit]=get_test_tray_RPDPM(y_test,feat_names,ages_test,timeline,...
    Sigma,A,B,C,D,G,loss_type, fit_type, optim_options,btstrp) 
%% Robust Parametric DPM Testing

% Estimating test subject parameters using biomaker values and ages
[K2, I2, J2] = size(y_test); % number of biomarkers, subjects, and visits
y_test_fit = zeros(K2, I2, J2, btstrp); % estimated values of biomarkers
y_tray_fit = zeros(K2, I2, size(timeline,2), btstrp); % estimated values of biomarkers

DPS_test = zeros(I2, J2, btstrp);
DPS_tray = zeros(I2, size(timeline,2), btstrp);

Alpha_test = zeros(I2, btstrp);
Beta_test = zeros(I2, btstrp);


for n = 1 : btstrp
    [Alpha_test(:, n), Beta_test(:, n), DPS_test(:, :, n),DPS_tray(:, :, n), y_test_fit(:, :, :, n),y_tray_fit(:,:,:,n)] =...
        rpdpm_test_tray(y_test, ages_test,timeline, ...
        Alpha_test(:, n), Beta_test(:, n), Sigma(:, n), A(:, n), B(:, n), C(:, n), D(:, n), G(:, n), loss_type, fit_type, optim_options);

%     [Alpha_test(:, n), Beta_test(:, n), DPS_test(:, :, n), y_test_fit(:, :, :, n)] = rpdpm_test(y_test, ages_test, ...
%         Alpha_test(:, n), Beta_test(:, n), Sigma(:, n), A(:, n), B(:, n), C(:, n), D(:, n), G(:, n), loss_type, fit_type, optim_options);
    
%     [Alpha_test(:, n), Beta_test(:, n), DPS_test(:, :, n), y_test_fit(:, :, :, n),...
%         tSample,tTrayectories] = ...
%         rpdpm_individual(y_test, ages_test, Alpha_test(:, n), Beta_test(:, n), ...
%         Sigma(:, n), A(:, n), B(:, n), C(:, n), D(:, n), G(:, n),...
%         loss_type, fit_type, optim_options,...
%         RID,feat_names);
end

% Test set modeling performance
NMAE_test = zeros(btstrp, 1); % overall test error per bootstrap
MAE_test = zeros(K2, btstrp); % test error per biomarker per bootstrap
for n = 1 : btstrp
    ae_test = abs(y_test - y_test_fit(:, :, :, n)); % absolute errors
    MAE_test(1 : K2, n) = nanmean(ae_test(1 : K2, :), 2); % mean absolute error
    ae_test = ae_test(1 : K2, :) ./ Sigma(1 : K2, n); % normalized (scaled) absolute errors
    % for k=1:K2
    %     fprintf('MAE (%s): %.3f\n',feat_names(k),mean(ae_test(k,:),'omitnan'));
    % end
    NMAE_test(n) = nanmean(ae_test(:)); % normalized mean absolute error
end

for k=1:K2
        fprintf('MAE (%s): %.3f\t',feat_names(k),mean(MAE_test(k,:,end),'omitnan'));
end

nmae_test_avg = mean(NMAE_test); % average of normalized mean absolute errors
nmae_test_std = std(NMAE_test); % standard deviation of normalized mean absolute errors
%mae_test_avg = mean(MAE_test, 2); % average of mean absolute errors
%mae_test_std = std(MAE_test, 0, 2); % standard deviation of mean absolute errors
fprintf('\nNMAE = %4.4f \x00B1 %4.4f \n', [nmae_test_avg, nmae_test_std]); % display test modeling performance

end
