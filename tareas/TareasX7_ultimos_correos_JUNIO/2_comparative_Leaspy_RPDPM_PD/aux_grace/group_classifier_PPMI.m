function acc = group_classifier_PPMI(time_onset,group,thresh_time)

labels=unique(group);
mask=group==labels(1);
time_NC=time_onset(mask,2);

mask=group==labels(2);

time_PD_bsl=time_onset(mask,1);
time_PD_end=time_onset(mask,2);


% acc=[sum(time_NC<=thresh_time)/length(time_NC),...
%      sum(time_PD>=thresh_time)/length(time_PD)];
% acc=[sum(time_NC<=thresh_time)/length(time_NC),...
%      sum(time_PD_bsl<=thresh_time & time_PD_end>=thresh_time)/length(time_PD_end)];
acc=[sum(time_NC<thresh_time)/length(time_NC),...
     sum(time_PD_end>=thresh_time)/length(time_PD_end)];

end