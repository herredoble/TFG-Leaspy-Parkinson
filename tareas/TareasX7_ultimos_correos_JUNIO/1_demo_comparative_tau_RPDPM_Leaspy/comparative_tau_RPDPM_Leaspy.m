
function comparative_tau_RPDPM_Leaspy

%% Input data file
addpath('./aux_grace');
addpath('./rpdpm'); % add path to the source codes
typeFeatures={'ROI','CAM','ROICAM','all'};
n_typeFeat= 2; %CAM
type_features=cellstr(typeFeatures{n_typeFeat});
bDeleteM3_M6_M18=false;
[ADNI,~,~,~,~] = getDataADNIMERGE_v1_2c(type_features,bDeleteM3_M6_M18);
ADNI_bsl=ADNI(ADNI.M==0,:);
%% Tau
tTau_Hyae=readtable('./data_multi/tau_Hyae.xlsx');
tTau_Hyae.ID=str2double(tTau_Hyae.ID);
load('./data_multi/tSample_tTray_AMC_tau_10_0y','tau');
tau(tau==inf)=nan;
tau_subj=mean(tau(:,2:end),2,'omitnan');
tau_subj=[tau(:,1),tau_subj];
%%
age_dem=nan(size(tau_subj,1),6);
for i=1:size(tau_subj,1)
    mask_subj=ADNI_bsl.RID==tau_subj(i,1);
    age_dem_gt = ADNI_bsl.AGE(mask_subj)-(ADNI_bsl.onset(mask_subj)/12);
    group=ADNI_bsl.stMCI(mask_subj);
    convert=ADNI_bsl.Convert(mask_subj);    
    try
        age_dem(i,:)=[tau_subj(i,1),group,convert,age_dem_gt,tau_subj(i,2),...
            tTau_Hyae.tau(tTau_Hyae.ID==tau_subj(i,1))];
    catch
        fprintf('RID %d without tau\n',tau_subj(i,1));
    end
end

%% Plot
figure(1);clf
mask_sMCI= age_dem(:,2)==1 & age_dem(:,3)==0;
mask_pMCI= age_dem(:,2)==1 & age_dem(:,3)==1;
mask_AD= age_dem(:,2)==2 & age_dem(:,3)==0;
subplot(1,2,1);
mask_nonan= ~isnan(age_dem(:,5));
plot(age_dem(mask_sMCI & mask_nonan,4),age_dem(mask_sMCI & mask_nonan,5),'g+');hold on;
plot(age_dem(mask_pMCI & mask_nonan,4),age_dem(mask_pMCI & mask_nonan,5),'b*');
plot(age_dem(mask_AD & mask_nonan,4),age_dem(mask_AD & mask_nonan,5),'ro');
legend('sMCI','pMCI','Dementia');
xlabel('Dementia onset [years]');
ylabel('Estimation [years]');
title('RPDPM');
corr=corrcoef(age_dem(mask_nonan,4),age_dem(mask_nonan,5));
fprintf('RPDPM dementia onsete correlation: %.3f\n',corr(1,2));

subplot(1,2,2);
mask_nonan= ~isnan(age_dem(:,6));
plot(age_dem(mask_sMCI & mask_nonan,4),age_dem(mask_sMCI & mask_nonan,6),'g+');hold on;
plot(age_dem(mask_pMCI & mask_nonan,4),age_dem(mask_pMCI & mask_nonan,6),'b*');
plot(age_dem(mask_AD & mask_nonan,4),age_dem(mask_AD & mask_nonan,6),'ro');
legend('sMCI','pMCI','Dementia');
xlabel('Dementia onset [years]');
ylabel('Estimation [years]');
title('Leaspy');
corr=corrcoef(age_dem(mask_nonan,4),age_dem(mask_nonan,6));
fprintf('Leaspy dementia onsete correlation: %.3f\n',corr(1,2));


end