function comparative_tau_RPDPM_Leaspy_PD
criterion_age =true;

% RPDPM MCI onset correlation: 0.976
% RPDPM MCI onset correlation with MCI/PDD subjects: 0.977
% RPDPM MCI onset correlation: 0.974
% RPDPM MCI onset correlation with MCI/PDD subjects: 0.978
% Leaspy MCI onset correlation: 0.028
% Leaspy MCI onset correlation with MCI/PDD subjects: 0.076
% 
% RPDPM MCI onset correlation: 0.495
% RPDPM MCI onset correlation with MCI/PDD subjects: 0.404
% RPDPM MCI onset correlation: 0.488
% RPDPM MCI onset correlation with MCI/PDD subjects: 0.679
% Leaspy MCI onset correlation: 0.322
% Leaspy MCI onset correlation with MCI/PDD subjects: 0.319

%% Input data file

addpath('./aux_grace');
n_dataset=1; %PPMIApril20
n_typeFeat=5; % typeFeatures={'NM','DAT','CSF','MOTOR','all'};
[PPMI,~,~,~,~] = getDataPPMI_MCI_v1_2(n_typeFeat,n_dataset);
PPMI_bsl = PPMI(PPMI.YEAR==0,:);


%% Tau
tTau_Leaspy_test=readtable('./data/IPs_media_test_tau.xlsx');
tTau_Leaspy_train=readtable('./data/IPs_media_train_tau.xlsx');
tTau_Leaspy=[tTau_Leaspy_test;tTau_Leaspy_train];
tTau_Leaspy.RID=str2double(tTau_Leaspy.ID);

% load('./data/tSample_tTray_MGSSCAU_tau_0y', 'tau');
% load('./data/tSample_tTray_JIS_tau_10_5y','tau');
load('./data/tSample_tTray_JIS_tau_10','tau');
tau(tau==inf)=nan;
tau_subj=mean(tau(:,2:end),2,'omitnan');
tau_subj=[tau(:,1),tau_subj];
%%
age_MCI=nan(size(tau_subj,1),6);
for i=1:size(tau_subj,1)
    mask_subj=PPMI_bsl.PATNO==tau_subj(i,1);
    age_MCI_gt = PPMI_bsl.age(mask_subj)-PPMI_bsl.onset(mask_subj);
    group=PPMI_bsl.cogstate_end(mask_subj);
    convert=PPMI_bsl.cogstate_end(mask_subj)>1;
    age_bsl = PPMI_bsl.age(mask_subj);
    try
        if(criterion_age == true)
            age_MCI(i,:)=[tau_subj(i,1),group,convert,age_MCI_gt,tau_subj(i,2),...
                tTau_Leaspy.tau(tTau_Leaspy.RID==tau_subj(i,1))];
        else
            age_MCI(i,:)=[tau_subj(i,1),group,convert,age_bsl-age_MCI_gt,age_bsl-tau_subj(i,2),...
                age_bsl-tTau_Leaspy.tau(tTau_Leaspy.RID==tau_subj(i,1))];
        end

    catch
        fprintf('RID %d without tau\n',tau_subj(i,1));
    end
end

%% Plot
figure(1);clf
mask_sCU= age_MCI(:,2)==1;
mask_pCU= age_MCI(:,2)==2;
mask_PDD= age_MCI(:,2)==3;
mask_MCI= age_MCI(:,2)>1;
subplot(1,2,1);
mask_nonan= ~isnan(age_MCI(:,5));
plot(age_MCI(mask_sCU & mask_nonan,4),age_MCI(mask_sCU & mask_nonan,5),'g+');hold on;
plot(age_MCI(mask_pCU & mask_nonan,4),age_MCI(mask_pCU & mask_nonan,5),'b*');
plot(age_MCI(mask_PDD & mask_nonan,4),age_MCI(mask_PDD & mask_nonan,5),'ro');
legend('CU','MCI','PDD');
xlabel('MCI onset [years]');
ylabel('Estimation [years]');
title('RPDPM');
corr=corrcoef(age_MCI(mask_nonan,4),age_MCI(mask_nonan,5));
fprintf('RPDPM MCI onset correlation: %.3f\n',corr(1,2));
corr=corrcoef(age_MCI(mask_MCI & mask_nonan,4),age_MCI(mask_MCI & mask_nonan,5));
fprintf('RPDPM MCI onset correlation with MCI/PDD subjects: %.3f\n',corr(1,2));


subplot(1,2,2);
mask_nonan= ~isnan(age_MCI(:,6));
plot(age_MCI(mask_sCU & mask_nonan,4),age_MCI(mask_sCU & mask_nonan,6),'g+');hold on;
plot(age_MCI(mask_pCU & mask_nonan,4),age_MCI(mask_pCU & mask_nonan,6),'b*');
plot(age_MCI(mask_PDD & mask_nonan,4),age_MCI(mask_PDD & mask_nonan,6),'ro');
legend('CU','MCI','PDD');
xlabel('MCI onset [years]');
ylabel('Estimation [years]');
title('Leaspy');
corr=corrcoef(age_MCI(mask_nonan,4),age_MCI(mask_nonan,6));
fprintf('Leaspy MCI onset correlation: %.3f\n',corr(1,2));

corr=corrcoef(age_MCI(mask_MCI & mask_nonan,4),age_MCI(mask_MCI & mask_nonan,6));
fprintf('Leaspy MCI onset correlation with MCI/PDD subjects: %.3f\n',corr(1,2));

end