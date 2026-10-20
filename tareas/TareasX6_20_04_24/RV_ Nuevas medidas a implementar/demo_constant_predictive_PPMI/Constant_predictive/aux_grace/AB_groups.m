function convertAB=AB_groups(grace)
convertAB=nan(size(grace,1),1);
[PD_mas,PD_menos]=Abeta_PD_PPMI;
RID=unique(PD_mas.PATNO);
for i=1:length(RID)
    mask=grace.RID==RID(i);
    if(sum(mask)>0)
        convertAB(mask)=1;
    else
        fprintf(2,'Error RID %d without in grace\n',RID(i));
    end
        
end

RID=unique(PD_menos.PATNO);
for i=1:length(RID)
    mask=grace.RID==RID(i);
    if(sum(mask)>0)
        convertAB(mask)=0;
    else
        fprintf(2,'Error RID %d without in grace\n',RID(i));
    end
end





end