function [PD_MCI,PD_HC]=cognitivePD_PPMI

%addpath('./selected_aux');

load('./data/PPMIApr20_5Years','PPMIApr20');

mask_PD=PPMIApr20.APPRDX==1;
PPMI_PD=PPMIApr20(mask_PD,:);
ID=unique(PPMI_PD.PATNO);
% PD_MCI_1=[];
% PD_MCI_2=[];
% PD_MCI_3=[];
mask_PD_MCI=false(size(PPMI_PD,1),1);
convertTime=zeros(size(PPMI_PD,1),1);
for i=1:length(ID)
    index=find(PPMI_PD.PATNO==ID(i));
    PD_subj=PPMI_PD(index,:);
    cogstate=PD_subj.cogstate; % 1(HC),2(MCI), 3(Dementia)
    MCI_testscores=PD_subj.MCI_testscores; % 0 (HC), 1(MCI)
    np1cog=PD_subj.NP1COG;  %0 (HC), 1(Slight), 2(Mild), 3 (Moderate), 4(Severate)
    EVENT_ID=PD_subj.EVENT_ID;
    data_ID=[vis2double(EVENT_ID),cogstate,MCI_testscores,np1cog];
    %% Criteria 1
    criteria1=0;
    cogstate=cogstate(~isnan(cogstate));
    if(sum(cogstate>1)>0 && cogstate(end)>1)
        %PD_MCI_1=[PD_MCI_1;ID(i)];
        criteria1=1;
    end
    %% Criteria 2
    criteria2=0;
    MCI_testscores=MCI_testscores(~isnan(MCI_testscores));
    if(sum(MCI_testscores==1)>0 && MCI_testscores(end)==1)
        %PD_MCI_2=[PD_MCI_2;ID(i)];
        criteria2=1;
    end
    %% Criteria 3
    criteria3=0;
    np1cog=np1cog(~isnan(np1cog));
    if(sum(np1cog>0)>1 && np1cog(end)>0)
        %PD_MCI_3=[PD_MCI_2;ID(i)];
        criteria3=1;
    end
    %% Criteria join
    if((criteria1+criteria3)>1)
%         PD_MCI=[PD_MCI;ID(i)];
%         %ID(i)
%         data_ID
%         pause;
        mask_PD_MCI(index)=true;
        idx_first_cr1=find(data_ID(:,2)>1,1);
        idx_first_cr3=find(data_ID(:,4)>0,1);
        convertTime(index)=(data_ID(idx_first_cr1,1)+data_ID(idx_first_cr3,1))/2;
    else
        convertTime(index)=data_ID(end,1);
        
    end
    
        
        
end

PPMI_PD.convertTime=convertTime;
Convert=zeros(size(PPMI_PD,1),1);
Convert(mask_PD_MCI)=1;
PPMI_PD.Convert=Convert;

PD_MCI=PPMI_PD(mask_PD_MCI,:);
PD_HC =PPMI_PD(~mask_PD_MCI,:);


%% delete
%rmpath('./selected_aux');

end





