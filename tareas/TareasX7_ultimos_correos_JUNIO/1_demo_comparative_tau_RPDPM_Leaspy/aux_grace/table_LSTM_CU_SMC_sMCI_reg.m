function [ADNI,matrix_LSTM,matrix_LSTM_wNaN,feature_names,percent,point_feat]=...
    table_LSTM_CU_SMC_sMCI_reg(batches,b_interp,tbl_delta_sCU_SMC_MCI,bOnly_pCU_pSMC)
close all
clc;


%addpath('./aux_grace');

%% data
n_typeFeat= 2; %CAM
typeFeatures={'ROI','CAM','ROICAM','all'};
type_features=cellstr(typeFeatures{n_typeFeat});
bDeleteM3_M6_M18=true;
[ADNI,~,~,~,~] = getDataADNIMERGE_v1_2c(type_features,bDeleteM3_M6_M18);

mask_CN=(ADNI.stMCI==0) | (ADNI.stMCI<0); % CU and SMC
mask_sMCI=ADNI.stMCI==1 & ADNI.Convert==0 & ADNI.DX_bl=='EMCI';
mask=mask_CN | mask_sMCI;

ADNI=ADNI(mask,:);
ADNI=sortrows(ADNI,{'RID','M'});
ADNI(ADNI.M==3,:)=[];
ADNI(ADNI.M==6,:)=[];
ADNI(ADNI.M==18,:)=[];


RID=unique(ADNI.RID);
n_lotes=max(2,batches);
if(b_interp)
    for i=1:length(RID)
        mask=ADNI.RID==RID(i);
        if(sum(mask)<n_lotes)
            ADNI(ADNI.RID==RID(i),:)=[];
        end
    end
end

DX_num=nan(size(ADNI,1),1);
DX_num(ADNI.DX=='CN')=0;
DX_num(ADNI.DX=='MCI')=1;
DX_num(ADNI.DX=='Dementia')=1;
ADNI.DX_num=DX_num;


feature_tbl_names=string(ADNI.Properties.VariableNames)';

selectedFeatures=[29;25;23;34]; % Learning- A13+ CDRSB+ TRABSCOR
change_sign_feat=false(length(selectedFeatures),1);
change_sign_feat(1)=true; %Learning-        

n_features=length(selectedFeatures);
%DX=61;
feature_names=feature_tbl_names(selectedFeatures);
features=ADNI(:,selectedFeatures).Variables;

%% stable subject selections with estimation time shift
[mask_cohort,ADNI]=add_delta_stable_subject(ADNI,tbl_delta_sCU_SMC_MCI,...
    bOnly_pCU_pSMC);
features=features(mask_cohort,:);

%% Incresead features
group=ADNI.stMCI>0;
for i=1:length(change_sign_feat)
    if(change_sign_feat(i)==true)
        features(:,i)=-features(:,i);
    end
end
% Transforming
[features,percent,point_feat]=percentil_interp(features,group);



%% LSTM matrix
RID=unique(ADNI.RID);
matrix_LSTM=[];
matrix_LSTM_wNaN=[];
for i=1:length(RID)
    fprintf('ID %d: %d\n',i,RID(i));
    tbl_subj=ADNI(ADNI.RID==RID(i),:);
    fea_subj=features(ADNI.RID==RID(i),:);
    feat_check=tbl_subj(:,selectedFeatures).Variables;
    years=tbl_subj.M/12;
    timeline=(0:max(years))';
    Convert=tbl_subj.Convert(1);
    if(RID(i)==926)
        fprintf('Error\n');       
    end    
%     if(RID(i)==4499)
%         fprintf('Error\n');
%         continue;
%     end
    if(length(timeline)~=length(years))       
        idx = find(~ismember(timeline, years));
        data_subj=nan(length(timeline),n_features);
        diag_subj=nan(length(timeline),1);
        onset_subj=nan(length(timeline),1);
        years_check=diag_subj;
        mask_samples=true(length(timeline),1);
        mask_samples(idx)=false;
        
        %data_subj(mask_samples,:)=tbl_subj(:,selectedFeatures).Variables;
        data_subj(mask_samples,:)=fea_subj;
        
        diag_subj_tmp =tbl_subj(:,end).Variables;
        if(Convert==1)
            idx_conv=find(diag_subj_tmp==1,1);
            diag_subj_tmp(idx_conv:end)=1;
        else
           group_bsl=diag_subj_tmp(1);
           diag_subj_tmp(:)=group_bsl;
        end
        diag_subj(mask_samples) =diag_subj_tmp;
        
        onset_subj(mask_samples,:)= tbl_subj.onset/12;
        
        years_check(mask_samples)=years;
    else
        %data_subj=tbl_subj(:,selectedFeatures).Variables;
        data_subj=fea_subj;
        diag_subj=tbl_subj(:,end).Variables;
        if(Convert==1)
            idx_conv=find(diag_subj==1,1);
            diag_subj(idx_conv:end)=1;
        else
           group_bsl=diag_subj(1);
           diag_subj(:)=group_bsl;
        end
            
        onset_subj = tbl_subj.onset/12;
    end
    
    %% Interpolation
    if(b_interp)
        data_subj_interp = get_Data_Interp_Subj(data_subj);
        onset_subj_interp = get_Onset_Interp_Subj(onset_subj); 
        diag_subj_interp = get_Diag_Interp_Subj(diag_subj,onset_subj_interp);
                      
    end
    
    
    
    
    
    %% Matrix n_bath x batch*n_features
    matrix_features_batch=get_flat_data_batch(data_subj,batches);
    if(size(matrix_features_batch,1)>1)
        diag_subj_temp=diag_subj(1:batches);
        diag_subj_temp=diag_subj_temp(~isnan(diag_subj_temp));
        diag_subj=[diag_subj_temp(end);diag_subj(batches+1:end)];
        
        onset_subj_temp=onset_subj(1:batches);
        onset_subj_temp=onset_subj_temp(~isnan(onset_subj_temp));
        onset_subj=[onset_subj_temp(end);onset_subj(batches+1:end)];
        
    else
        diag_subj=diag_subj(~isnan(diag_subj));
        diag_subj=diag_subj(end);
        
        onset_subj=onset_subj(~isnan(onset_subj));
        onset_subj=onset_subj(end);
       
        
    end
    
    if(b_interp)
        matrix_features_batch_wNaN=get_flat_data_batch(data_subj_interp,batches);
        diag_subj_interp=diag_subj_interp(batches:end);
        onset_subj_interp=onset_subj_interp(batches:end);
    end
        
    
   
    n_batch=size(matrix_features_batch,1);
    matrix_features_batch=[repmat(RID(i),n_batch,1),matrix_features_batch,...
        diag_subj,repmat(Convert,n_batch,1),timeline(batches:end),onset_subj];
    
    if(size(matrix_features_batch,2)==(5+(n_features*batches)))
           matrix_LSTM=[matrix_LSTM;matrix_features_batch];
    else
        fprinf('Error\n');
    end
    
    if(b_interp)
        matrix_features_batch_wNaN=[repmat(RID(i),n_batch,1),matrix_features_batch_wNaN,...
            diag_subj_interp,repmat(Convert,n_batch,1),timeline(batches:end),onset_subj_interp];
        if(size(matrix_features_batch_wNaN,2)==(5+(n_features*batches)))
            matrix_LSTM_wNaN=[matrix_LSTM_wNaN;matrix_features_batch_wNaN];
        else
            fprinf('Error\n');
        end
            
    end
        
                
       


end

        



end


function [mask_cohort,ADNI_out]=add_delta_stable_subject(ADNI,tbl_delta_sCU_SMC_MCI,bOnly_pCU_pSMC)
mask_cohort=false(size(ADNI,1),1);
mask_cohort(ADNI.Convert==1)=true;

if(bOnly_pCU_pSMC==false)
    RID=tbl_delta_sCU_SMC_MCI.RID;
    for i=1:length(RID)
        mask_subj=ADNI.RID==RID(i);
        if(sum(mask_subj)>0)
            tmp=ADNI.onset(mask_subj);
            ADNI.onset(mask_subj)=ADNI.onset(mask_subj)+(tbl_delta_sCU_SMC_MCI.delta(i)*12);
            tmp2=ADNI.onset(mask_subj);
            %checking stable subjects
            group=ADNI.stMCI(mask_subj)>=1;
            if(group(1)~=tbl_delta_sCU_SMC_MCI.Group(i))
                fprintf('Error delta estimation in RID %d\n',...
                    tbl_delta_sCU_SMC_MCI.RID(i));
            end
            if((group(1)==1) && (tbl_delta_sCU_SMC_MCI.delta(i)>0))
                mask_cohort(mask_subj)=true;
            elseif((group(1)==0) && (tmp2(end)<0))
                mask_cohort(mask_subj)=true;
            else
            end
         end
    end
end

ADNI_out=ADNI(mask_cohort,:);
end



function onset_subj_interp = get_Onset_Interp_Subj(onset_subj)

mask_onset=~isnan(onset_subj);

if(sum(mask_onset)~=length(onset_subj))
    timeline=(0:length(onset_subj)-1)';
    onset_subj_interp=onset_subj;
    try
        onset_interp = interp1(timeline(mask_onset),onset_subj(mask_onset),...
            timeline(~mask_onset),'linear','extrap');%, 'spline');
    catch
        onset_interp=zeros(sum(~mask_onset),1);
        %onset_interp(:)=onset_subj(mask_onset);
    end
    onset_subj_interp(~mask_onset)=onset_interp;
else
   onset_subj_interp= onset_subj;
end



end

function data_subj_interp = get_Data_Interp_Subj(data_subj)

n_features=size(data_subj,2);
timeline=0:size(data_subj,1)-1;
data_subj_interp=data_subj;
for j=1:n_features
    feature=data_subj(:,j);
    mask_feature=~isnan(feature);
    if(sum(mask_feature)~=length(feature))
        try
            feat_interp = interp1(timeline(mask_feature),feature(mask_feature),...
                timeline(~mask_feature),'linear','extrap');%, 'spline');
        catch
            feat_interp=zeros(sum(~mask_feature),1);
            feat_interp(:)=feature(mask_feature);
        end
        feature(~mask_feature)=feat_interp;
        data_subj_interp(:,j)=feature;
    end
     
end

end

function diag_subj_interp = get_Diag_Interp_Subj(diag_subj, onset_interp )

diag_subj_interp=diag_subj;
mask_diag=~isnan(diag_subj);
if(sum(mask_diag)~=length(diag_subj))
    labels=unique(diag_subj);
    if(length(labels)==1)
        diag_subj_interp(~mask_diag)=labels;
    else
        idx=find(mask_diag==0);
        diag_subj_interp(idx)=onset_interp(idx)>0;
%         for i=1:length(idx)
%             if(idx(i)>1)
%                 diag_subj_interp(idx(i))=diag_subj_interp(idx(i)-1);
%             else
%                 diag_subj_interp(idx(i))=labels(1);
%             end            
%         end
        
    end
        
    
end


end

function flattened_data=get_flat_data_batch(data,batches)
points=size(data,1);
feature_dim=size(data,2);

if(points<=batches)
    points=1;
    flattened_data = nan(points, batches * feature_dim);
    flattened_data(1, 1:size(data, 1) * feature_dim) = reshape(data', 1, []);   
else
    points=points-batches+1;
    flattened_data = nan(points, batches * feature_dim);
    rows=1:batches;
    for i = 1:points
        d_tmp=data(rows,:)';
        flattened_data(i,:) = reshape(d_tmp(:), 1, []);
        rows=rows+1;
    end

    
end
    
end

        