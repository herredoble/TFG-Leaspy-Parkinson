function [mask_NC,Converts,convertTime,onset]=convertTime_NC_v1(ADNIMERGE)

% function [mask_NC,Converts,convertTime,onset,Converts_amiloide,convertTime_amiloide,onset_amiloide]=...
%     convertTime_NC_amiloide(ADNIMERGE)

 
mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');
mask_NC=(ADNIMERGE.DX_bl=='CN') & mask_DX; %NC
ADNIMERGE=ADNIMERGE(mask_NC,:); %CU:2952v 515s, NC: 2771v 409s, (old)4233v 512s
%RID=unique(ADNIMERGE.RID);

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
% % ADNI_pNC= ADNIMERGE(Converts_amiloide>0,:); %pNC: 867v 121s
% % ADNI_sNC= ADNIMERGE(Converts_amiloide==0,:); %sNC: 1492v 215s

%% Clinical converts
mask_pNC= (ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI'); % people who convert to MCI/AD
Converts=false(size(ADNIMERGE,1),1);
pNC_subjs=unique(ADNIMERGE.RID(mask_pNC));
for i=1:length(pNC_subjs)
    index=find(ADNIMERGE.RID==pNC_subjs(i));
    DX=ADNIMERGE.DX(index);
    if((DX(end)=='Dementia') || (DX(end)=='MCI'))
        Converts(index)=true;
    end 
end

% ADNIMERGE.Convert=Converts;
% ADNI_pNC= ADNIMERGE(Converts,:); %pNC: 721v 81s
% ADNI_sNC= ADNIMERGE(Converts==0,:); %sNC: 1638v 255s

%% pCU (amiloide)
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



%% pCU (clinical) 810v 96s
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

%% sCU (amiloide)
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



%% sCU (clinical) 419s 2142v
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


