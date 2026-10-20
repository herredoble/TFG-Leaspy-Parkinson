function grace=set_tableGrace_PPMI(PPMI,features,raw,show)



%% Table for R
RID=unique(PPMI.PATNO);

%% Table for R
num_features=size(features,2);

grace=zeros(length(features(:)),7);

index_begin=1;
for i=1:length(RID)
    index=find(PPMI.ID==RID(i));
    index_end=index_begin+(length(index)*num_features)-1;
    features_subj=features(index,:);
    raw_subj=raw(index,:);
    id_marker=zeros(size(features_subj));
    for j=1:num_features
        id_marker(:,j)=j;
    end
    [t,index_time]=sort(PPMI.t(index));
    features_subj=features_subj(index_time,:);
    raw_subj=raw_subj(index_time,:);
    grace(index_begin:index_end,1)=repmat(t/12,num_features,1);    
    grace(index_begin:index_end,2)=features_subj(:);
    grace(index_begin:index_end,3)=raw_subj(:);   
    grace(index_begin:index_end,4)=id_marker(:);
    grace(index_begin:index_end,5)=repmat(PPMI.ID(index),num_features,1);
    grace(index_begin:index_end,6)=repmat(PPMI.Convert(index),num_features,1);
    grace(index_begin:index_end,7)=repmat(PPMI.PATNO(index),num_features,1);
    grace(index_begin:index_end,8)=repmat(PPMI.convertTime(index),num_features,1);
    
    index_begin=index_end+1;
    
end

grace=grace(~isnan(grace(:,2)),:);
grace=array2table(grace,'VariableNames',...
    {'t','Y','raw','Outcome','id','Convert','RID','ConvertTime'});

%% Correlations
if(show)
    outcomes=unique(grace.Outcome);
    LMEmodels=cell(length(outcomes),1);
    grace=cleanID_sinVisit(grace);

    for i=1:length(outcomes)
        grace_outcome=grace(grace.Outcome==outcomes(i),:);
        LMEmodels{i}=getLME_outcome(grace_outcome);
    end
    [inter,slope]=corrOutcomeLME(LMEmodels);

    % corr_inter=corrcov(cov(inter));
    % corr_slope=corrcov(cov(slope));
    
    figure(size(features,2)+1);clf;
    %subplot(1,2,1);
    corrplot(array2table(inter,'VariableNames',feature_names));
    %title('MCI test scores (Intercept)');
    title('Cogstate (Intercept)');
   
    figure(size(features,2)+2);clf;
    %subplot(1,2,2);
    corrplot(array2table(slope,'VariableNames',feature_names));
    %title('MCI test scores (Slope)');
    title('Cogstate (Slope)');
end    


end



