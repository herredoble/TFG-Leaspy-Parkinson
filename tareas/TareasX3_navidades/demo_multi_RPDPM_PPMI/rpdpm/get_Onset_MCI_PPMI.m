function [ages_onset,gt_onset,metrics]=get_Onset_MCI_PPMI(DPS_test,posterior_bayes,ages_test,labels_test,btstrp)

I=size(DPS_test,1);
ages_onset=nan(I,btstrp);
for n = 1 : btstrp
    dps_test_n = DPS_test(:, :, n);
    prob_MCI=posterior_bayes{2, n}(dps_test_n)>.5;
    for i=1: I
        if(sum(prob_MCI(i,:)>0))
            mask_values= ~isnan(ages_test(i,:));
            age_subj=ages_test(i,mask_values);
            label_subj=prob_MCI(i,mask_values);
            idx=find(label_subj);
            if(label_subj(end))
                if(idx(1)>1)
                    ages_onset(i,n)=(age_subj(idx(1))+age_subj(idx(1)-1))/2;
                else
                    ages_onset(i,n)=age_subj(idx(1));
                end
            end
        end %convert
    end %subjects
end %btstrp

gt_onset=get_GT_onset(labels_test,ages_test);

%% metrics
% PD_CU
metrics=table;
percentage_PD_CU=nan(btstrp,1);
percentage_PD_MCI=nan(btstrp,1);
corr_PD_MCI_age=nan(btstrp,1);
corr_PD_MCI_reserve=nan(btstrp,1);

numPD_CU=sum(isnan(gt_onset));

for n = 1 : btstrp
    percentage_PD_CU(n)=sum(isnan(ages_onset(:,n)) & isnan(gt_onset))/numPD_CU*100;
end
metrics.percentage_PD_CU=percentage_PD_CU;

% PD_MCI_PDD
age_bsl=ages_test(:,1);
numPD_MCI=sum(~isnan(gt_onset));

for n = 1 : btstrp
    percentage_PD_MCI(n)=sum(~isnan(ages_onset(:,n)) & ~isnan(gt_onset))/numPD_MCI*100;
    cov=corrcoef(rmmissing([ages_onset(:,n),gt_onset]));
    corr_PD_MCI_age(n)=cov(1,2);
    cov=corrcoef(rmmissing([ages_onset(:,n)-age_bsl,gt_onset-age_bsl]));
    corr_PD_MCI_reserve(n)=cov(1,2);
    
end
metrics.percentage_PD_MCI=percentage_PD_MCI;
metrics.corr_PD_MCI_age=corr_PD_MCI_age;
metrics.corr_PD_MCI_reserve=corr_PD_MCI_reserve;
%% show
fprintf('Percentage PD-CU: %.1f (%.2f)\n',mean(metrics.percentage_PD_CU),...
    std(metrics.percentage_PD_CU));
fprintf('Percentage PD-MCI/PDD: %.1f (%.2f)\n',mean(metrics.percentage_PD_MCI),...
    std(metrics.percentage_PD_MCI));
fprintf('Correlation age conversion: %.2f (%.2f)\n',mean(metrics.corr_PD_MCI_age),...
    std(metrics.corr_PD_MCI_age));
fprintf('Correlation cognition reserve: %.2f (%.2f)\n',mean(metrics.corr_PD_MCI_reserve),...
    std(metrics.corr_PD_MCI_reserve));

end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gt_onset=get_GT_onset(labels_test,ages_test)
I=size(labels_test,1);
gt_onset=nan(I,1);
MCI=strcmp(labels_test,'MCI');

for i=1: I
    if(sum(MCI(i,:)>0))
        mask_values= ~isnan(ages_test(i,:));
        age_subj=ages_test(i,mask_values);
        label_subj=MCI(i,mask_values);
        idx=find(label_subj);
        if(label_subj(end))
            if(idx(1)>1)
                gt_onset(i)=(age_subj(idx(1))+age_subj(idx(1)-1))/2;
            else
                gt_onset(i)=age_subj(idx(1));
            end
        end
    end %convert
end %subjects

end