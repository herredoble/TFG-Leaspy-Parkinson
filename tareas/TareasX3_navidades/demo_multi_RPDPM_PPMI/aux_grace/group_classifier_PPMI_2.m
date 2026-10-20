function acc = group_classifier_PPMI_2(time_onset,group,thresh_time)

% mask=group=='PD-MCI-';
% time_sNC=time_onset(mask,2);
% 
% mask=group=='PD-MCI+';

mask=group=='PD-CU';
time_sNC=time_onset(mask,2);

mask= (group=='PD-MCI') | (group=='PDD');

time_pNC=[time_onset(mask,1),time_onset(mask,2)];

acc=[sum(time_sNC<thresh_time)/length(time_sNC),...
     sum(time_pNC(:,2)>thresh_time)/size(time_pNC,1),...
     sum(time_pNC(:,1)<thresh_time)/size(time_pNC,1)];

end