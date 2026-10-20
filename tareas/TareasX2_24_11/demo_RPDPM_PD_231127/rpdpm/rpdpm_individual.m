function [alpha, beta, s, y_fit,tSample,tTrayectories] = rpdpm_individual(y, t, alpha0, beta0,...
    sigma, a, b, c, d, g, loss_type, fit_type, optims,...
    RID,name_biomarkers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 Source Code Author Information                  %%%%%
%%%%%                     Mostafa Mehdipour Ghazi                     %%%%%
%%%%%                   mostafa.mehdipour@gmail.com                   %%%%%
%%%%%                      Created on 01/08/2017                      %%%%%
%%%%%                      Updated on 01/08/2020                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  This code is an implementation of the algorithm published in:  %%%%%
%%%%%  Robust parametric modeling of Alzheimer's disease progression  %%%%%
%%%%%                https://arxiv.org/abs/1908.05338                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing the robust parametric disease progression model

% I: number of subjects
% K: number of biomarkers
% J: number of visits (time points)
% y: biomarker measurements [K x I x J]
% t: visiting age of subjects [I x J]
% [alpha0, beta0]: initial values of subject-specific parameters each of which [I x 1]
% sigma: scale factor of biomarker errors [K x 1]
% [a, b, c, d, g]: biomarker-specific parameters each of which [K x 1]
% loss_type: robust estimation loss function type
% fit_type: fitting function type
% optims: optimization options for fitting biomarkers
% [alpha, beta]: subject-specific parameters each of which [I x 1]
% s: disease progression scores (DPS) [I x J]
% y_fit: estimated biomarker measurements [K x I x J]

% Initializing the model parameters
[K, I, J] = size(y);
alpha = zeros(I, 1);
beta = zeros(I, 1);

% Assigning ranges of the model parameters
lb = [zeros(I, 1); - Inf(I, 1)]; % lower bounds of parameters [alpha, beta]
ub = [Inf(I, 1); Inf(I, 1)]; % upper bounds of parameters [alpha, beta]

% Weights for normalizing the residuals
w = zeros(I, J, K);
for i = 1 : I
    w(i, :) = 1 / sum(sum(~isnan(y(:, i, :)), 3), 1); % normalized w.r.t. the number of available points per subject
end
w_i = permute(w, [1, 3, 2]);

% Replicating biomarker-specific parameters for all time-points
a_i = repmat(a, 1, J);
b_i = repmat(b, 1, J);
c_i = repmat(c, 1, J);
d_i = repmat(d, 1, J);
g_i = repmat(g, 1, J);
sigma_i = repmat(sigma, 1, J);

% Estimating subject-specific parameters
for i = 1 : I
    
    % Finding measures with available values
    y_i = squeeze(y(:, i, :));
    t_i = repmat(t(i, :), K, 1);
    idx = find(~isnan(t_i(:)) & ~isnan(y_i(:)));
    
    % Fitting sigmoids to biomarkers
    x0 = [alpha0(i); beta0(i)];
    F = @(x) objective_subject(x, y_i(idx), t_i(idx), w_i(i, idx)', sigma_i(idx), a_i(idx), b_i(idx), c_i(idx), d_i(idx), g_i(idx), loss_type, fit_type);
    params = lsqnonlin(F, x0, lb(i : I : end), ub(i : I : end), optims);
    
    % Updating subject-specific parameters
    alpha(i) = params(1);
    beta(i) = params(2);
    
end

% Estimating DPS values and measurements

% s = repmat(beta, 1, J) + repmat(alpha, 1, J) .* t; % linear DPS
% y_fit = NaN(K, I, J);
% for k = 1 : K
%     for i = 1 : I
%         y_fit(k, i, :) = sigmoid_function(fit_type, s(i, :), a(k), b(k), c(k), d(k), g(k));
%     end
% end

% s = repmat(beta, 1, J) + repmat(alpha, 1, J) .* t; % linear DPS
% s_t = repmat(beta(i), 1, J2) + repmat(alpha(i), 1, J2) .* timeline; % linear DPS
% y_fit = NaN(K, I, J);
% y_timeline = NaN(K, I, J2);
% 
% for k = 1 : K
%     for i = 1 : I
%         y_fit(k, i, :) = sigmoid_function(fit_type, s(i, :), a(k), b(k), c(k), d(k), g(k));
%         y_timeline(k, i, :) = sigmoid_function(fit_type, s_t, a(k), b(k), c(k), d(k), g(k));        
%     end
% end

s = repmat(beta, 1, J) + repmat(alpha, 1, J) .* t; % linear DPS
y_fit = NaN(K, I, J);

tSample=table;
tTrayectories=table;
for i = 1 : I
    clear tSubj;
    tSubj=table;
    mask_dps= ~isnan(s(i,:));
    tSubj.RID=repmat(RID(i),sum(mask_dps),1);
    tSubj.DPS=s(i,mask_dps)';
    tSubj.Age=t(i,mask_dps)';
    
    timeline=-5+min(t(i,:)):0.5:5+max(t(i,:));
    J2=length(timeline);
    s_t = repmat(beta(i), 1, J2) + repmat(alpha(i), 1, J2) .* timeline; % linear DPS
    clear tSubTray;
    tSubTray=table;
    tSubTray.RID=repmat(RID(i),J2,1);
    tSubTray.DPS=s_t';
    tSubTray.Age=timeline';

    for k = 1 : K
        y_fit(k, i, :) = sigmoid_function(fit_type, s(i, :), a(k), b(k), c(k), d(k), g(k));
        f = sigmoid_function(fit_type, s_t, a(k), b(k), c(k), d(k), g(k));
        ySample=squeeze(y(k, i, :))';
        yfit=squeeze(y_fit(k, i, :))';
        tSubj=[tSubj,table(ySample(mask_dps)','VariableNames',name_biomarkers(k))];
        tSubj=[tSubj,table(yfit(mask_dps)','VariableNames',string(sprintf('est-%s',name_biomarkers(k))))];
        tSubTray=[tSubTray,table(f','VariableNames',string(sprintf('tray-%s',name_biomarkers(k))))];
    end
    
    tSample=[tSample;tSubj];
    tTrayectories=[tTrayectories;tSubTray];
    
end

%plot_tayectory_individual(alpha(1),beta(1),fit_type, a, b, c, d, g, squeeze(y(:,1,:)),squeeze(t(1,:)))
end



% function plot_tayectory_individual(alpha,beta,fit_type, a, b, c, d, g,y,t)
% J=length(t);
% K=size(y,1);
% timeline=-5+min(t):0.5:5+max(t);
% J2=length(timeline);
% mask_samples= ~isnan(t);
% s = repmat(beta, 1, J2) + repmat(alpha, 1, J2) .* timeline; % linear DPS
% t_ = repmat(beta, 1, J) + repmat(alpha, 1, J) .* t; % linear DPS
% 
% for k = 1 : K
%     figure;
%     f = sigmoid_function(fit_type, s, a(k), b(k), c(k), d(k), g(k)); % fit logistic curve
%     plot(t_(mask_samples),y(k,mask_samples),'+');hold on;
%     plot(s,f,'*');hold on;
% end
% 
% 
% end