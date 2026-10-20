function [mask_MCI,Converts,convertTime,onset]= convertTime_MCI_v1(ADNIMERGE)

% function [mask_MCI,Converts,convertTime,onset,Converts_amiloide,convertTime_amiloide,onset_amiloide]=...
%     convertTime_MCI_amiloide(ADNIMERGE)

 

mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');

% mask_MCI= ADNIMERGE.DX_bl=='LMCI' & mask_DX;
% mask_MCI= (ADNIMERGE.DX_bl=='EMCI') & mask_DX;
mask_MCI=((ADNIMERGE.DX_bl=='LMCI') | (ADNIMERGE.DX_bl=='EMCI')) & mask_DX;

%LMCI: 3474v 542s: 1920v 261s_MCI 1554v 281s_pMCI
%eMCI: v s:  v 262s_MCI v 46s_pMCI
%elMCI: v s: v 523s_MCI v 327s_pMCI
ADNIMERGE=ADNIMERGE(mask_MCI,:); 
% RID=unique(ADNIMERGE.RID); % 1015s 5694v 

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

%% Clinical converts
mask_pMCI= (ADNIMERGE.DX=='Dementia'); % people who convert to AD
Converts=false(size(ADNIMERGE,1),1);
pMCI_subjs=unique(ADNIMERGE.RID(mask_pMCI));
for i=1:length(pMCI_subjs)
    index=find(ADNIMERGE.RID==pMCI_subjs(i));
    DX=ADNIMERGE.DX(index);
    if(DX(end)=='Dementia')
        Converts(index)=true;
    end
end

% ADNI_pMCI= ADNIMERGE(Converts>0,:); %pMCI: 1768v 256s
% ADNI_sMCI= ADNIMERGE(Converts==0,:); %sMCI: 2731v 444s

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
% % ADNI_MCI= ADNIMERGE(Converts_amiloide>0,:); %MCI: 2579v 409s
% % ADNI_nMCI= ADNIMERGE(Converts_amiloide==0,:); %nMCI: 1920v 291s
% 


%% pMCI 345s 2304v
RID=unique(ADNIMERGE.RID(Converts==1));
convertTime=zeros(size(ADNIMERGE,1),1);
onset=convertTime;
for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
    viscode=ADNIMERGE.M(index);
    [viscode,ord]=sort(viscode);
    if(length(ord)>1)
        if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
            fprintf(2,'No sorting visits from RID %d.\n',RID(i));
        end
    end
    DX=ADNIMERGE.DX(index(ord));
    ct=find(DX=='Dementia',1);
    convertTime(index)=viscode(ct);
    time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
    if(ct>1)
        Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
            ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;
    else
        Delta_days=0;      
    end
    onset(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30;
    

end


%% MCI (amiloide)
% RID=unique(ADNIMERGE.RID(Converts_amiloide>0));
% convertTime_amiloide=nan(size(ADNIMERGE,1),1);
% onset_amiloide=convertTime_amiloide;
% for i=1:length(RID)
%     index=find(ADNIMERGE.RID==RID(i));
%     viscode=ADNIMERGE.M(index);
%     [viscode,ord]=sort(viscode);
%     if(length(ord)>1)
%         if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
%             fprintf(2,'No sorting visits from RID %d.\n',RID(i));
%         end
%     end
%     if(Converts_amiloide(index(1))==1)
%         AV45=ADNIMERGE.AV45(index(ord));
%         ct=find(AV45>=AV45_th,1);
%         
%         convertTime_amiloide(index)=viscode(ct);
%         time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
%         if(ct>1)          
%             Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
%                 ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;            
%         else
%             Delta_days=0;
%         end
%         onset_amiloide(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30;
%        
%     elseif(Converts_amiloide(index(1))==2)
%         pTauAB=ADNIMERGE.PTAU(index(ord))./ADNIMERGE.ABETA(index(ord));
%         ct=find(pTauAB>=pTauAB_th,1);
%         
%         convertTime_amiloide(index)=viscode(ct);
%         time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
%         if(ct>1)          
%             Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
%                 ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;            
%         else
%             Delta_days=0;
%         end
%         
%         onset_amiloide(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30;
%        
%     else
%     end
%     
% end


%% sMCI
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

%% nMCI (amiloide)
% RID=unique(ADNIMERGE.RID(Converts_amiloide==0));
% for i=1:length(RID)
%     index=find(ADNIMERGE.RID==RID(i));
%     viscode=ADNIMERGE.M(index);
%     [viscode,ord]=sort(viscode);
%     if(length(ord)>1)
%         if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
%             fprintf(2,'No sorting visits from RID %d.\n',RID(i));
%         end
%     end    
%     convertTime_amiloide(index)=viscode(end);
%     onset_amiloide(index)=-viscode(end);
% end






end


