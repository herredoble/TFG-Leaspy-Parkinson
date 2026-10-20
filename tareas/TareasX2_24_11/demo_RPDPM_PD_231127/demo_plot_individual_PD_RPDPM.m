function demo_plot_individual_PD_RPDPM

close all;
clc;

load('./data/tau_tray_sample_10_LSHSMA','tau','tSample','tTrayec');
load('./data/percentil_LSHSMA.mat', 'percent', 'point_feat');

RID=unique(tSample.RID);
markers_name=["ABETA","HVLTRDLY","LNS","MEANCAU","SDMTOTAL","SFT"];
change_sign_feat=true(length(markers_name),1);


bSalir=false;
idx=0;
while(bSalir==false)
    reply=str2double(input('Idx from patient (1:423) or 0 (exit):','s'));
    if(isnan(reply))
        idx=idx+1;
        fprintf('%d\n',idx);
    else
        idx=reply;
    end
    if(idx~=0)
        plot_individual_original_LSHSMA(tau(idx,:),tSample(tSample.RID==RID(idx),:),...
            tTrayec(tTrayec.RID==RID(idx),:),markers_name);
        plot_individual_norm_LSHSMA(tau(idx,:),tSample(tSample.RID==RID(idx),:),...
            tTrayec(tTrayec.RID==RID(idx),:),markers_name,...
            change_sign_feat,percent,point_feat);
        
    else
       bSalir=true; 
    end
end

end


function plot_individual_original_LSHSMA(tau,tSample,tTrayec,markers_name)
RID=tSample.RID(1);

%t_onset
meanAge_MCI_RPDPM=mean(tau(~isinf(tau)));
% stdAge_MCI_RPDPM=mean(tau(~isinf(tau)));
% if(isnan(meanAge_MCI_RPDPM))
%     meanAge_MCI_RPDPM=max(tTrayec.Age);
%     stdAge_MCI_RPDPM=0;
% end

%Clinical Diagnose
DX=string(tSample.DX);
DX=DX(~strcmp(DX,""));
diag=string(unique(DX));
if(numel(diag)==2)
    diag='pCU';
end
mask_CU=strcmp(tSample.DX,'CN');
mask_MCI=strcmp(tSample.DX,'MCI');
mask_NaN= (mask_CU==0) & (mask_MCI==0);
    
%Markers    
tS_col_names=string(tSample.Properties.VariableNames)';
tS_col_bstrp=extractBefore(tS_col_names,'_');
tT_col_names=string(tTrayec.Properties.VariableNames)';
tT_col_bstrp=extractBefore(tT_col_names,'_');

figure(1);clf;
for k=1:numel(markers_name)
    col_S=find(tS_col_names==markers_name(k),1);
    cols_SB=find(tS_col_bstrp==markers_name(k));
    cols_TB=find(tT_col_bstrp==markers_name(k));
    subplot(2,3,k);
    samples=table2array(tSample(:,col_S));
    estimate=table2array(tSample(:,cols_SB));
    mask_samples = ~isnan(samples);
    if(sum(mask_samples)>0)
        mae_bstrp=calculateMAE(samples(mask_samples),estimate(mask_samples,:));
    else
        mae_bstrp=nan;
    end
    plot(tSample.Age(mask_CU & mask_samples),samples(mask_CU & mask_samples),'b*');hold on;
    plot(tSample.Age(mask_MCI & mask_samples),samples(mask_MCI & mask_samples),'ro');
    plot(tSample.Age(mask_NaN & mask_samples),samples(mask_NaN & mask_samples),'kp');
    max_v=max(samples);
    min_v=min(samples);
    trayect=table2array(tTrayec(:,cols_TB));
    for n=1:length(cols_TB)
        plot(tTrayec.Age,trayect(:,n),'-');
    end
    max_v=max([max_v;trayect(:)]);
    min_v=min([min_v;trayect(:)]);      
    
    min_v=min_v-(0.2*min_v);
    max_v=max_v+(0.2*max_v);
    marker_th=linspace(min_v,max_v,10)';
    %t_onset
    plot(repmat(meanAge_MCI_RPDPM,length(marker_th),1),marker_th,'--',...
        'LineWidth',2);
    %plot(repmat(meanAge_MCI_RPDPM+(3*stdAge_MCI_RPDPM),length(marker_th),1),marker_th,'--');
    %plot(repmat(meanAge_MCI_RPDPM-(3*stdAge_MCI_RPDPM),length(marker_th),1),marker_th,'--');

    ylim([min_v,max_v]);
    
    grid on;
    if(sum(~isnan(mae_bstrp))>0)
        title(sprintf('RID: %d DX: %s MAE = %.2f(%.2f)',RID,diag,...
            mean(mae_bstrp),std(mae_bstrp)));
    else 
        title(sprintf('RID: %d DX: %s)',RID,diag));
    end

    xlabel('years')
    ylabel(markers_name(k));
        
    if(strcmp(diag,'CN'))
        legend('CU','Estimate');
    elseif(strcmp(diag,'MCI'))
        legend('MCI','Estimate');
    else
        legend('CU','MCI','Estimate');
    end
end


end

function mae_bstrp=calculateMAE(samples,estimate)
bstrp=size(estimate,2);
mae_bstrp=mean(abs(repmat(samples,1,bstrp)-estimate))';
end

function plot_individual_norm_LSHSMA(tau,tSample,tTrayec,markers_name,...
    change_sign_feat,percent,point_feat)
% function plot_individual_norm_LACT(tau,tSample,tTrayec,markers_name,...
%     change_sign_feat,percent,point_feat,...
%      A,B,C,D,G)
RID=unique(tTrayec.RID);
% fit_type = 'proposed'; % fitting function type

%t_onset
meanAge_MCI_RPDPM=mean(tau(~isinf(tau)));
%stdAge_MCI_RPDPM=mean(tau(~isinf(tau)));
if(isnan(meanAge_MCI_RPDPM))
    meanAge_MCI_RPDPM=max(tTrayec.Age);
    %stdAge_MCI_RPDPM=0;
end

%Clinical Diagnose
DX=string(tSample.DX);
DX=DX(~strcmp(DX,""));
diag=string(unique(DX));
if(numel(diag)==2)
    diag='pCU';
end
    
%Markers    
tT_col_names=string(tTrayec.Properties.VariableNames)';
tT_col_bstrp=extractBefore(tT_col_names,'_');

figure(2);clf;
for k=1:numel(markers_name)
    cols_TB=find(tT_col_bstrp==markers_name(k));
    trayect=table2array(tTrayec(:,cols_TB));
    trayect=mean(trayect,2,'omitnan');
    if(change_sign_feat(k)==false)
        norm_trayect=interp1(point_feat{k},percent{k},trayect);
    else
        norm_trayect=interp1(point_feat{k},percent{k},-trayect);
    end
    
    plot(tTrayec.Age,norm_trayect,'-'); hold on;
    %plot(tTrayec.dps_1,norm_trayect,'-');
%     btstrp=10;
%     f_avg=zeros(size(tTrayec,1),1);
%     for n = 1 : btstrp
%         f = sigmoid_function(fit_type, tTrayec.dps_2, A(k, n), B(k, n), C(k, n), D(k, n), G(k, n)); % fit logistic curve
%         f = (f - D(k, n)) / (A(k, n) - D(k, n)); % scale to [0 1]
%         f_avg = f_avg + f;
%     end
%     f_avg=f_avg/btstrp;
    %plot(tTrayec.Age,f_avg,':','LineWidth',2);
    
end

marker_th=0:0.1:1;
plot(repmat(meanAge_MCI_RPDPM,length(marker_th),1),marker_th,'--',...
        'LineWidth',2);
ylim([0,1]);
grid on;
title(sprintf('RID: %d DX: %s',RID,diag));
xlabel('years')
legend(markers_name);       
end

