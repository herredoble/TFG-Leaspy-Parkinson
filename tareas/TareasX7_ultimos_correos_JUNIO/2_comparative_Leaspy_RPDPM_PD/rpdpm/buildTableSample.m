function tSample=buildTableSample(RID,y,y_fit,ages,DPS,name_biomarkers)


K=length(name_biomarkers);
bstrp=size(DPS,3);
tSample=table;
for i = 1 : length(RID)
    clear tSubj;
    tSubj=table;
    mask_dps= ~isnan(ages(i,:));
    tSubj.RID=repmat(RID(i),sum(mask_dps),1);
    tSubj.Age=ages(i,mask_dps)';
    
    if(~isempty(y))
        for k = 1 : K
            ySample=squeeze(y(k, i, :))';
            tSubj=[tSubj,table(ySample(mask_dps)','VariableNames',name_biomarkers(k))];
        end
    end
        
    
    for n=1:bstrp
        time_DSP=squeeze(DPS(i, :,n));
        tSubj=[tSubj,table(time_DSP(mask_dps)','VariableNames',string(sprintf('dps_%d',n)))];
        
        for k = 1 : K
            y_k_i_n=squeeze(y_fit(k, i, :,n))';
            tSubj=[tSubj,table(y_k_i_n(mask_dps)','VariableNames',...
                string(sprintf('%s_%d',name_biomarkers(k),n)))];
        end
        
        
    end
    
    tSample=[tSample;tSubj];
   
end



end