function [mask_AD,Converts,convertTime,onset]= convertTime_AD_v1(ADNIMERGE)

% function [mask_AD,Converts,convertTime,onset,Converts_amiloide,convertTime_amiloide,onset_amiloide]= ...
%     convertTime_AD_amiloide(ADNIMERGE)

 

mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');

mask_AD=(ADNIMERGE.DX_bl=='AD') & mask_DX;

%AD: 3474v 542s: 1920v 261s_MCI 1554v 281s_pMCI
ADNI_AD=ADNIMERGE(mask_AD,:); 
%RID=unique(ADNIMERGE.RID); %395s 1192v
%% sort
% RID=unique(ADNIMERGE.RID);
% index_all=[];
% for i=1:length(RID)
%     index=find(ADNIMERGE.RID==RID(i));
%     [val,idx]=sort(ADNIMERGE.Month(index));
%     if((val(1)~=0) || (ADNIMERGE.VISCODE(index(idx(1)))~='bl'))
%         fprintf('Error %d\n',RID(i));
%     else
%         index_all=[index_all;index(idx)];
%     end
%     
% end
% ADNIMERGE=ADNIMERGE(index_all,:);

%% clinical regresion
Converts=false(size(ADNIMERGE,1),1);

mask_regresion= ADNI_AD.DX=='MCI' | ADNI_AD.DX=='CN';
rSubjs=unique(ADNI_AD.RID(mask_regresion));
for i=1:length(rSubjs)
    index=find(ADNI_AD.RID==rSubjs(i));
    DX=ADNI_AD.DX(index);
    if(DX(end)=='MCI' || DX(end)=='CN')
        index=find(ADNIMERGE.RID==rSubjs(i));
        Converts(index)=true;
    end
end


mask_AD(Converts)= 0;
convertTime=zeros(sum(mask_AD),1);
Converts=convertTime;
onset=convertTime;

%% sDementia
ADNIMERGE=ADNIMERGE(mask_AD,:); 
RID=unique(ADNIMERGE.RID(Converts==0));

for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
    viscode=ADNIMERGE.M(index);
    [viscode,ord]=sort(viscode);
    if(length(ord)>1)
        if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
            fprintf(2,'No sorting visits from RID %d.\n',RID(i));
        end
    end    
    convertTime(index)=viscode(end);
    onset(index)=-viscode(end);

end

% ADNI_rAD= ADNIMERGE(Converts>0,:); %AD: 8v 2s
% ADNI_AD= ADNIMERGE(Converts==0,:); %AD: 752v 231s


%% Amiloide converts
% RID=unique(ADNIMERGE.RID);
% Converts_amiloide=zeros(size(ADNIMERGE,1),1);
% 
% AV45_th=1.16;
% % TauAB_th=0.33;
% pTauAB_th=0.028;
% 
% for i=1:length(RID)
%     index=find(ADNIMERGE.RID==RID(i));
%     AV45=ADNIMERGE.AV45(index);
%     AV45=AV45(~isnan(AV45));    
% %     TauAB=ADNIMERGE.TAU(index)./ADNIMERGE.ABETA(index);
% %     TauAB=TauAB(~isnan(TauAB));
%     pTauAB=ADNIMERGE.PTAU(index)./ADNIMERGE.ABETA(index);
%     pTauAB=pTauAB(~isnan(pTauAB));
%     
%     if(~isempty(AV45))
%         if(AV45(end)>=AV45_th)
%             Converts_amiloide(index)=1;
%         end
%     elseif(~isempty(pTauAB))
%         if(pTauAB(end)>=pTauAB_th)
%             Converts_amiloide(index)=2;
%         end
%     else
%         Converts_amiloide(index)=nan;
%     end
%     
% end
% 
% % ADNI_AD= ADNIMERGE(Converts_amiloide>0,:); %AD: 655v 198s
% % ADNI_nAD= ADNIMERGE(Converts_amiloide==0,:); %nAD: 105v 35s
% 
% 
% 
% convertTime_amiloide=convertTime;
% onset_amiloide=convertTime;
% 



end


