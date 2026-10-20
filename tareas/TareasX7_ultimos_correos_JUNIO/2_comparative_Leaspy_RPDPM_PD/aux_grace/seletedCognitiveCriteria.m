function [converts,convertTime]=seletedCognitiveCriteria(grace,n_cogCriteria)

converts=nan(size(grace,1),1);
convertTime=nan(size(grace,1),1);
[PD_MCImas,PD_MCImenos]=cognitivePD_PPMI;
PPMI=[PD_MCImas;PD_MCImenos];
[PPMI,ID_sinCog]=Cognitive_groups(PPMI);
%% Cognitive criteria
cognitiveCriteria={'MCI_test','UPDRS1','SITE'};
if(strcmp('MCI_test',cognitiveCriteria{n_cogCriteria}))
    converts_=PPMI.MCI_testscores_end;
    convertTime_=PPMI.convertTimeMCI_test;
elseif(strcmp('UPDRS1',cognitiveCriteria{n_cogCriteria}))
    converts_=PPMI.np1cog_end>1; %0 (HC), 1(Slight), 2(Mild), 3 (Moderate), 4(Severate)
    convertTime_=PPMI.convertTimeNp1cog;    
elseif(strcmp('SITE',cognitiveCriteria{n_cogCriteria}))
    PPMI=estimateCogState(PPMI,ID_sinCog);
    %converts_=double(PPMI.cogstate_end>1); %1(HC),2(MCI), 3(Dementia)
    converts_=PPMI.cogstate_end; %1(HC),2(MCI), 3(Dementia)
    mask_nan=isnan(PPMI.cogstate_end);
    converts_(mask_nan)=nan;
    convertTime_=PPMI.convertTimeCogstate;
    
else
end
    


RID=PPMI.PATNO(PPMI.YEAR==0);


for i=1:length(RID)
    mask_grace=grace.RID==RID(i);
    index=find(PPMI.PATNO==RID(i));
    
    converts(mask_grace)=converts_(index(1));
    convertTime(mask_grace)=convertTime_(index(1));
        
end



end