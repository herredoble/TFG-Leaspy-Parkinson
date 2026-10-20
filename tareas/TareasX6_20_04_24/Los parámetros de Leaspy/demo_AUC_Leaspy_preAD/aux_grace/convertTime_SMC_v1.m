function [mask_SMC,Converts,convertTime,onset]= convertTime_SMC_v1(ADNIMERGE)

% function [mask_SMC,Converts,convertTime,onset,Converts_amiloide,convertTime_amiloide,onset_amiloide]=...
%     convertTime_SMC_v1(ADNIMERGE)

 
mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');
mask_SMC=(ADNIMERGE.DX_bl=='SMC') & mask_DX; %NC
ADNIMERGE=ADNIMERGE(mask_SMC,:); %SMC: 456v 105s
% RID=unique(ADNIMERGE.RID); 759V 303s

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

%% Amiloide converts
% Converts_amiloide=zeros(size(ADNIMERGE,1),1);
% 
% AV45_th=1.11;
% % TauAB_th=0.33;
% % pTauAB_th=0.028;
% pTauAB_th=-880;
% 
% for i=1:length(RID)
%     index=find(ADNIMERGE.RID==RID(i));
%     AV45=ADNIMERGE.AV45(index);
%     AV45=AV45(~isnan(AV45));    
% %     TauAB=ADNIMERGE.TAU(index)./ADNIMERGE.ABETA(index);
% %     TauAB=TauAB(~isnan(TauAB));
% 
% %     pTauAB=ADNIMERGE.PTAU(index)./ADNIMERGE.ABETA(index);
%     pTauAB=-ADNIMERGE.ABETA(index);
%     pTauAB=pTauAB(~isnan(pTauAB));
%     
% %     if(~isempty(AV45))
% %         if(AV45(end)>=AV45_th)
% %             Converts_amiloide(index)=1;
% %         end
% %     elseif(~isempty(pTauAB))
% %         if(pTauAB(end)>=pTauAB_th)
% %             Converts_amiloide(index)=2;
% %         end
% %     else
% %         Converts_amiloide(index)=nan;
% %     end
% 
%     AV45plus=false;
%     if(~isempty(AV45))
%         AV45plus=max(AV45)>=AV45_th;
%         if(AV45plus==true)
%             Converts_amiloide(index)=1;
%         end        
%     end
%     if(~isempty(pTauAB))
%         pTauplus=max(pTauAB)>=pTauAB_th;       
%         if((pTauplus==true) && (AV45plus==true))
%             Converts_amiloide(index)=3;
%         elseif(pTauplus==true)
%             Converts_amiloide(index)=2;
%         end
%     end
%     if(isempty(AV45) && isempty(pTauAB))
%         Converts_amiloide(index)=nan;
%     end
% end


%% Clinical converts
mask_pSMC= (ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI'); % people who convert to MCI/AD
Converts=false(size(ADNIMERGE,1),1);
pSMC_subjs=unique(ADNIMERGE.RID(mask_pSMC));
for i=1:length(pSMC_subjs)
    index=find(ADNIMERGE.RID==pSMC_subjs(i));
    DX=ADNIMERGE.DX(index);
    if((DX(end)=='Dementia') || (DX(end)=='MCI'))
        Converts(index)=true;
    end 
end


%% pSMC (amiloide)
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
% 


%% pSMC (clinical) 22s 91 v
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
    ct=find((DX=='MCI') | (DX=='Dementia'),1);
    convertTime(index)=viscode(ct);
    time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
    Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
        ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;
    onset(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30;
    
end

%% sSMC (amiloide) 281s 668v
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



%% sSMC (clinical)
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





end


