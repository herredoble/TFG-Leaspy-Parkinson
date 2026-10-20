function summaryPPMI_bsl(raw_bsl,group,feature_names,PPMI)

gr=unique(group);
num_gr=length(gr);

% if(num_gr>3)
%     num_gr=3;
%     gr=gr(1:3);
% end

if(num_gr==2)
    p_test=true;
    %p_test=false;
else
    p_test=false;
end

% fprintf('Subjects: ');
% for i=1:num_gr
%     fprintf(' %d (%s), ',sum(group==gr(i)),char(gr(i)));
% end

fprintf('Subjects: ');
for i=1:num_gr
    mask_gr=group==gr(i);
    N=sum(mask_gr);
    Female=sum(PPMI.gen(mask_gr)==2);%Female
    mask_APOE=~isnan(PPMI.APOE_e4) & mask_gr;
    APOE4=sum(PPMI.APOE_e4(mask_APOE)>0);
    fprintf('%s: N=%d , Age= %.1f (%.1f), (%.1f %.1f)\t',char(gr(i)),N,mean(PPMI.age(mask_gr)),std(PPMI.age(mask_gr)),...
        min(PPMI.age(mask_gr)),max(PPMI.age(mask_gr)));
    fprintf('Female: %d (%.1f)\n',Female,Female/N*100);
    fprintf('APOE4: %d  %d (%.1f)\n',APOE4,sum(mask_APOE),APOE4/sum(mask_APOE)*100);
    fprintf('Education: %.2f (%.2f)\n',mean(PPMI.EDUCYRS(mask_gr)),...
        std(PPMI.EDUCYRS(mask_gr)));
end

fprintf('\n\n');


for i=1:numel(feature_names)
    fprintf('%s: ',feature_names{i});
    for j=1:num_gr
        fprintf('%s (%d): %.2f (%.2f) (%.2f %.2f),',char(gr(j)),sum(isnan(raw_bsl(group==gr(j),i))==0),...
            mean(raw_bsl(group==gr(j),i),'omitnan'),std(raw_bsl(group==gr(j),i),'omitnan'),...
            min(raw_bsl(group==gr(j),i)),max(raw_bsl(group==gr(j),i)));
    end
    if(p_test)
        try
            fprintf(' pvalue: %.4f',ranksum(raw_bsl(group==gr(1),i),raw_bsl(group==gr(2),i)));
        catch
            fprintf('No evaluable');
        end
    else
        p_cutoff=0.01;
        multcomparation(raw_bsl(:,i),group,p_cutoff);

    end
    fprintf('\n');

end

end

function multcomparation(raw,group,p_cutoff)
[p,~,stats] = anova1(raw(:,1),group,'off');
fprintf('& %.4f\t',p);
[results,~,~,gnames] = multcompare(stats,'CType','bonferroni','Display','off');
res=results(results(:,end)<p_cutoff,[1 2]);
gnames=string(gnames);
fprintf('& \\multirow{2}{*}{');
for i=1:size(res,1)
    fprintf('%s-%s,\t',gnames(res(i,1)),gnames(res(i,2)));
end
fprintf('}');
end

