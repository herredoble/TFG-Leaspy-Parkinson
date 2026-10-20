function [PPMI_NC,PPMI_PD]=selectedSubjectPPMI(dataset)


if(strcmp(dataset,'FS_507'))
    load('./data/Tresults_507v_175s_170820','Tresults');
elseif(strcmp(dataset,'PA_633'))
    load('./data/Tresults_633v_215s_31082020','Tresults');
end
Tresults.Month=Tresults.years*12;
PPMIMERGE=Tresults;
clear Tresults       



%% sort
RID=unique(PPMIMERGE.fsidbase);
index_all=[];
for i=1:length(RID)
    index=find(PPMIMERGE.fsidbase==RID(i));
    if(strcmp(dataset,'FS_507') || strcmp(dataset,'PA_633'))
        [val,idx]=sort(PPMIMERGE.Month(index));
        if((val(1)~=0) || (PPMIMERGE.VISCODE(index(idx(1)))~='bl'))
            fprintf('Error %d\n',RID(i));
        else
            index_all=[index_all;index(idx)];
        end
    end
    
end
PPMIMERGE=PPMIMERGE(index_all,:);


%% Converts
PPMIMERGE.Converts=strcmp(PPMIMERGE.diagnose,'PD');

PPMI_PD= PPMIMERGE(PPMIMERGE.Converts,:); 
PPMI_NC= PPMIMERGE(PPMIMERGE.Converts==0,:);



end