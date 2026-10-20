function Apos = getApos_PPMI(PPMI)

% abeta_cutoff=700;
abeta_cutoff=683.45;

Apos=nan(size(PPMI,1),1);
ID=unique(PPMI.PATNO);
for i=1:length(ID)
    mask_subj=PPMI.PATNO==ID(i);
    abeta=PPMI.abeta(mask_subj);
    if(~isnan(abeta(1)))
       abeta_bsl=abeta(1);
    elseif(sum(~isnan(abeta))>0)
        abeta=abeta(~isnan(abeta));
        abeta_bsl=abeta(1);
    else
        abeta_bsl=[];
    end
    if(~isempty(abeta_bsl))
        if(abeta_bsl<=abeta_cutoff)
            Apos(mask_subj)=1;
        else
            Apos(mask_subj)=0;
        end
    end
        
       
        
end




end