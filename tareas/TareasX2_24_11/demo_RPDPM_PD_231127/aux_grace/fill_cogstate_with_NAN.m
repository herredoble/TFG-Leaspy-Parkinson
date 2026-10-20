function [cogstate_interp,ID_sinCog] = fill_cogstate_with_NAN(PPMI)

ID_sinCog=subjects_sinCogstate(PPMI.PATNO,PPMI.cogstate);
cogstate=estimate_cogstate_sin(PPMI.PATNO,PPMI.cogstate,...
    PPMI.MCI_testscores, PPMI.NP1COG,ID_sinCog);

PPMI.cogstate=cogstate;
cogstate_interp=PPMI.cogstate;

ID=unique(PPMI.PATNO);
ID_sinCog=[];
for i=1:length(ID)
%     if(ID(i)==3007)
%         fprintf('Error\n');
%     end
    mask_subj=PPMI.PATNO==ID(i);
    tbl_subj=PPMI(mask_subj,:);
    cogstate_subj=tbl_subj.cogstate; % 1(HC),2(MCI), 3(Dementia)
    mask_cog= ~isnan(cogstate_subj);
    if(size(tbl_subj,1)==sum(mask_cog))
        continue;
    elseif(sum(mask_cog)==0)
        ID_sinCog=[ID_sinCog;ID(i)];
    else
        NP1COG_subj=tbl_subj.NP1COG;
        cogstat_interp_sub=interp_cogstate(cogstate_subj,NP1COG_subj);
        cogstate_interp(mask_subj)=cogstat_interp_sub;   
    end
        

end

end


function cogstate_est=estimate_cogstate_sin(PATNO,cogstate,...
    MCI_testscores,NP1COG,ID_sinCog)

cogstate_est=cogstate;
for i=1:length(ID_sinCog)
    mask_subj=PATNO==ID_sinCog(i);
    MCI_test_subj=MCI_testscores(mask_subj)+1;
    NP1COG_subj= (NP1COG(mask_subj)>1)+1;
    
    cog_subj=nan(sum(mask_subj),1);
    mask_MCI_test= (~isnan(MCI_test_subj));
    mask_NP1COG= (~isnan(NP1COG_subj));
    
    mask_equal_criteria= MCI_test_subj==NP1COG_subj & mask_MCI_test & mask_NP1COG;
    if(sum(mask_equal_criteria)>0)
        cog_subj(mask_equal_criteria)=MCI_test_subj(mask_equal_criteria);       
    else
        MCI_= MCI_test_subj(mask_MCI_test);
        NPI_= NP1COG_subj(mask_NP1COG);
        idx=find(mask_NP1COG);
        cog_subj(idx(1))=min(MCI_(1),NPI_(1));
        cog_subj(idx(end))=min(MCI_(end),NPI_(end));        
    end
    cogstate_est(mask_subj)=cog_subj;

end

end


function ID_sinCog=subjects_sinCogstate(PATNO,cogstate)
ID=unique(PATNO);
ID_sinCog=nan(length(ID),1);
for i=1:length(ID)
    cogstate_subj=cogstate(PATNO==ID(i));
    mask_cog= ~isnan(cogstate_subj);
    if(sum(mask_cog)==0)
        ID_sinCog(i)=ID(i);
    end   
end
ID_sinCog=ID_sinCog(~isnan(ID_sinCog));

end

function cogstat_interp=interp_cogstate(cogstate_subj,NP1COG_subj)
cogstat_interp=cogstate_subj;
mask_nan= ~isnan(cogstate_subj);
cogstate_subj_know=cogstate_subj(mask_nan);
diag=unique(cogstate_subj_know);
if(length(diag)==1)
    if(diag==1)
        cogstat_interp(mask_nan==0)=diag;
    else
        cogstat_interp(mask_nan==0)= (NP1COG_subj(mask_nan==0)>1) + 1;
        cogstat_interp(end)=cogstate_subj_know(end);
    end
elseif(cogstate_subj_know(end)==1)
    cogstat_interp(mask_nan==0)=1;
else
    idx=find(cogstat_interp>1,1);
    mask_aux=true(length(cogstat_interp),1);
    mask_aux(idx:end)=false;
    cogstat_interp(mask_aux & mask_nan==0)=1;
    cogstat_interp(mask_aux==0 & mask_nan==0)=cogstat_interp(idx);
    cogstat_interp(end)=cogstate_subj_know(end);
end
        
end



