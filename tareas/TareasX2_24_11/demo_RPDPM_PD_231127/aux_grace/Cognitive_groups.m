function [PPMI,ID_sinCog]=Cognitive_groups(PPMI)

ID=unique(PPMI.PATNO);
cogstate=nan(size(PPMI,1),1);
MCI_testscores=nan(size(PPMI,1),1);
np1cog=nan(size(PPMI,1),1);

convertTimeCogstate=nan(size(PPMI,1),1);
convertTimeMCI_test=nan(size(PPMI,1),1);
convertTimeNp1cog=nan(size(PPMI,1),1);

ID_sinCog=[];

for i=1:length(ID)
    index=find(PPMI.PATNO==ID(i));
    if(ID(i)==3027)
        fprintf('Error\n');
    end
    PD_subj=PPMI(index,:);
    cogstate_subj=PD_subj.cogstate; % 1(HC),2(MCI), 3(Dementia)
    MCI_testscores_subj=PD_subj.MCI_testscores; % 0 (HC), 1(MCI)
    np1cog_subj=PD_subj.NP1COG;  %0 (HC), 1(Slight), 2(Mild), 3 (Moderate), 4(Severate)
    year_subj=PD_subj.YEAR;
    %% Criteria cogstate: 1(HC),2(MCI), 3(Dementia)
    year_cog=year_subj(~isnan(cogstate_subj));
    cogstate_subj=cogstate_subj(~isnan(cogstate_subj));
    if(isempty(cogstate_subj)==0)
        cogstate(index)=cogstate_subj(end);
        if(cogstate_subj(end) > 1)
            convertTimeCogstate(index)=repmat(min(year_cog(cogstate_subj>1)),...
                length(index),1);
        else
            convertTimeCogstate(index)=repmat(max(year_cog),length(index),1);          
        end
    else
        %fprintf(2,'Subject %d without cogstate data\n',ID(i));
        ID_sinCog=[ID_sinCog;1 ID(i)];
    end
    %% Criteria MCI_testscores:  0 (HC), 1(MCI)
    year_MCI=year_subj(~isnan(MCI_testscores_subj));
    MCI_testscores_subj=MCI_testscores_subj(~isnan(MCI_testscores_subj));
    if(isempty(MCI_testscores_subj)==0)        
        MCI_testscores(index)=MCI_testscores_subj(end);
        if(MCI_testscores_subj(end) > 0)         
            convertTimeMCI_test(index)=repmat(min(year_MCI(MCI_testscores_subj>0)),...
                length(index),1);
        else
            convertTimeMCI_test(index)=repmat(max(year_MCI),length(index),1);          
        end
    else
        fprintf(2,'Subject %d without MCI testscores data\n',ID(i));
        ID_sinCog=[ID_sinCog;2 ID(i)];
    end
    %% Criteria np1cog: 0 (HC), 1(Slight), 2(Mild), 3 (Moderate), 4(Severate)
    year_cog=year_subj(~isnan(np1cog_subj));
    np1cog_subj=np1cog_subj(~isnan(np1cog_subj));
    if(isempty(np1cog_subj)==0)
        np1cog(index)=np1cog_subj(end);
        if(np1cog_subj(end) > 1)         
            convertTimeNp1cog(index)=repmat(min(year_cog(np1cog_subj>1)),...
                length(index),1);
        else
            convertTimeNp1cog(index)=repmat(max(year_cog),length(index),1);          
        end
    else
        fprintf(2,'Subject %d without np1cog data\n',ID(i));
        ID_sinCog=[ID_sinCog;3 ID(i)];
    end
end
    
PPMI.cogstate_end=cogstate;
PPMI.MCI_testscores_end=MCI_testscores;
PPMI.np1cog_end=np1cog;

PPMI.convertTimeCogstate=convertTimeCogstate;
PPMI.convertTimeMCI_test=convertTimeMCI_test;
PPMI.convertTimeNp1cog=convertTimeNp1cog;
%% delete

end





