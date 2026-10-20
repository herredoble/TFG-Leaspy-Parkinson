function [ages_onset,gt_onset]=get_Onset_MCI(DPS_test,posterior_bayes,ages_test,labels_test,btstrp)

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
            if(label_subj(idx(end)))
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
numPD_CU=sum(isnan(gt_onset));
percentage_PD_CU=nan(btstrp,1);
for n = 1 : btstrp
    percentage_PD_CU(n)=sum(isnan(ages_onset(:,n)) & isnan(gt_onset))/numPD_CU;
end

% PD_MCI_PDD
numPD_MCI=sum(~isnan(gt_onset));
percentage_PD_MCI=nan(btstrp,1);
for n = 1 : btstrp
    percentage_PD_MCI(n)=sum(~isnan(ages_onset(:,n)) & ~isnan(gt_onset))/numPD_MCI;
end





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
        if(label_subj(idx(end)))
            if(idx(1)>1)
                gt_onset(i)=(age_subj(idx(1))+age_subj(idx(1)-1))/2;
            else
                gt_onset(i)=age_subj(idx(1));
            end
        end
    end %convert
end %subjects

end