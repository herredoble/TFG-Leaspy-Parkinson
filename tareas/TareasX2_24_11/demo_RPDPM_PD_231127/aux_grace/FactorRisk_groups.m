function factorRisk_plus=FactorRisk_groups(grace)
factorRisk_plus=nan(size(grace,1),1);
[PD_MCImas,PD_MCImenos]=cognitivePD_PPMI;
PPMI=[PD_MCImas;PD_MCImenos];

%% Factor risk
upsit_bsl=PPMI.upsit(PPMI.YEAR==0);
% converts=PPMI.Convert(PPMI.YEAR==0);
% Mdl = fitcdiscr(upsit_bsl,converts);
thres_upsit=20; %66.7 MCI-, 61.5 MCI+
RID=PPMI.PATNO(PPMI.YEAR==0);


for i=1:length(RID)
    mask=grace.RID==RID(i);
    if(sum(mask)>0)
        if(upsit_bsl(i)<thres_upsit)
            factorRisk_plus(mask)=1;
        else
            factorRisk_plus(mask)=0;
        end
            
    else
        fprintf(2,'Error RID %d without in grace\n',RID(i));
    end
        
end


end