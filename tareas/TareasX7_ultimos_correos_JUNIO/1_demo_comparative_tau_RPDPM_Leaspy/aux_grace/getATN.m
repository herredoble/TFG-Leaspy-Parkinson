function ATN=getATN(ADNI)

% AB_th=880;
pTau_th = 27;
pTauAB_th=0.028;

ATN=nan(size(ADNI,1),1);


[ABETA_min,TAU_max,pTAU_max,pTauAB_max,FDG_min]=getCSF_ATN(ADNI);

ADNI.ABETA_min=ABETA_min;
ADNI.TAU_max=TAU_max;
ADNI.PTAU_max=pTAU_max;
ADNI.pTauAB_max=pTauAB_max;
ADNI.FDG_min=FDG_min;

% ATN(ADNI.ABETA_min>AB_th)=0; % A-
ATN(ADNI.pTauAB_max<pTauAB_th)=0; % A-

% ATN(ADNI.ABETA_min<=AB_th & ADNI.PTAU_max<pTau_th)=1; %A+T-
% ATN(ADNI.ABETA_min<=AB_th & ADNI.PTAU_max>=pTau_th)=2; %A+T+
ATN(ADNI.pTauAB_max>=pTauAB_th & ADNI.PTAU_max<pTau_th)=1; %A+T-
ATN(ADNI.pTauAB_max>=pTauAB_th & ADNI.PTAU_max>=pTau_th)=2; %A+T+




end


function [ABETA_min,TAU_max,pTAU_max,pTauAB_max,FDG_min]=getCSF_ATN(ADNIMERGE)
%% Amiloide converts
RID=unique(ADNIMERGE.RID);
ABETA_min=nan(size(ADNIMERGE,1),1);
TAU_max=nan(size(ADNIMERGE,1),1);
pTAU_max=nan(size(ADNIMERGE,1),1);
pTauAB_max=nan(size(ADNIMERGE,1),1);
FDG_min=nan(size(ADNIMERGE,1),1);


for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
    ABETA=ADNIMERGE.ABETA(index);
    ABETA=ABETA(~isnan(ABETA));    
    Tau=ADNIMERGE.TAU(index);
    Tau=Tau(~isnan(Tau));   
    pTau=ADNIMERGE.PTAU(index);
    pTau=pTau(~isnan(pTau));
    pTauAB=ADNIMERGE.PTAU(index)./ADNIMERGE.ABETA(index);
    pTauAB=pTauAB(~isnan(pTau));
    
    FDG=ADNIMERGE.FDG(index);
    FDG=FDG(~isnan(FDG));
    
    if(~isempty(ABETA))
        ABETA_min(index)=min(ABETA);
    end
    if(~isempty(Tau))
         TAU_max(index)=max(Tau);
    end
    if(~isempty(pTau))
         pTAU_max(index)=max(pTau);
    end
    if(~isempty(FDG))
        FDG_min(index)=min(FDG);
    end
    if(~isempty(pTauAB))
         pTauAB_max(index)=max(pTauAB);
    end
    
end

end

