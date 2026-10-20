function crosstab_analysis(group, feat_categorical, feature_name)

for i=1:size(feat_categorical,2)    
    [tb,~,p_val]=crosstab(group,feat_categorical(:,i));
    fprintf('%s: %.3f (p-value)\n',feature_name(i),p_val);
    n_cohort=sum(tb,2);
    prb =[tb(:,1)./n_cohort*100,tb(:,2)./n_cohort*100];
    lbs_group=string(unique(group));
    for j=1:length(lbs_group)
        fprintf('%s %s %d: %d (%.1f) %d (%.1f)\n',lbs_group(j),feature_name(i),...
            n_cohort(j),tb(j,1),prb(j,1),tb(j,2),prb(j,2));
            
    end
    
    [p_val,lbs]=mult_crosstab(group,feat_categorical(:,i));
    for j=1:length(p_val)
        fprintf('%s %s: %.4f\n',lbs(j,1),lbs(j,2),p_val(j));
    end
    
%     Tbl=table(group,feat_categorical(:,i),'VariableNames',["Group",feature_name(i)]);
%     heatmap(Tbl,feature_name(i),'Group');
%     pause
end

end

function [p_val,lbs]=mult_crosstab(group,var_categorical)
labels=unique(group);
n=length(labels);
n_comb= factorial(n)/2/factorial(n-2);
p_val=zeros(n_comb,1);
lbs=strings(n_comb,2);
idx=1;
for i=1:n-1
    for j=i+1:n
        mask=(group==labels(i)) | (group==labels(j));
        gr=categorical(group(mask)~=labels(i),[0,1],string([labels(i),labels(j)]));
        cat=var_categorical(mask);       
        [~,~,p,grp]=crosstab(gr,cat);
        p_val(idx)=p;
        lbs(idx,:)=string(grp(:,1))';
        idx=idx+1;
       
    end
end

end
