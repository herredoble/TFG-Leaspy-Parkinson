function [m_corr,m_scores] = measures_grace_PD_MCI(grace,dd1,percent,point_feat,...
    change_sign_feat,feature_names,original_scale,show)

 show_extra=false;

%% Number of visits
outcomes=unique(dd1.outcome);
num_max_visits=0;
for i=1:length(outcomes)
    num_visits=sum(dd1.outcome==outcomes(i));
    if(num_visits>num_max_visits)
        num_max_visits=num_visits;      
    end
end


%% MAE & t_app
MAE=nan(length(unique(grace.id))*2,length(outcomes));
t_app=nan(length(unique(grace.id)),3,length(outcomes));
yhat_app=nan(length(unique(grace.id)),3,length(outcomes));

grace_t_ghat= nan(num_max_visits,5,length(outcomes));



for i=1:length(outcomes)
    mask_outcome=dd1.outcome==outcomes(i);
    dd_outcome=dd1(mask_outcome,:);

    id=unique(dd_outcome.id);

    idx_first=1;
      
    for j=1:length(id)
        % Grace
        index=find(dd_outcome.id==id(j));
        idx=idx_first:(idx_first+length(index)-1);
        t=dd_outcome.argvals(index)+dd_outcome.gamma0(index);
        [t_max,i_max]=max(t);
        [t_min,i_min]=min(t);
        t_grace=[t_min,t_max];
        %yihat=dd_outcome.alpha0+(dd_outcome.alpha1.*dd_outcome.argvals)+dd_outcome.ghat;
        yihat=[dd_outcome.yihat(index(i_min)),dd_outcome.yihat(index(i_max))];
        y_grace=[dd_outcome.y(index(i_min)),dd_outcome.y(index(i_max))];
        yihat_grace=yihat;
 

        ghat_grace=dd_outcome.ghat(index);
        t_ghat_grace=t;

%         per2feat=interp1(percent{i},point_feat{i},yihat')';
%         per2feat_ghat_grace=interp1(percent{i},point_feat{i},ghat_grace')';
%         if(change_sign_feat(i))
%             per2feat=-per2feat;
%             per2feat_ghat_grace=-per2feat_ghat_grace;
%         end
%         yihat_grace_raw=per2feat;

        if(original_scale==true)
            per2feat=interp1(percent{i},point_feat{i},yihat')';
            per2feat_ghat_grace=interp1(percent{i},point_feat{i},ghat_grace')';
            if(change_sign_feat(i))
                per2feat=-per2feat;
                per2feat_ghat_grace=-per2feat_ghat_grace;
            end
            yihat_grace_raw=per2feat;
        else
            yihat_grace_raw=yihat;
            per2feat_ghat_grace=ghat_grace;           
        end

        
        grace_t_ghat(idx,1,i)=repmat(id(j),length(index),1);
        grace_t_ghat(idx,2,i)=t_ghat_grace;
        grace_t_ghat(idx,3,i)=per2feat_ghat_grace;
        
        idx_first=idx_first+length(index);

        
        
        MAE([2*j-1,2*j],i)=abs(y_grace-yihat_grace)';
        
        
        
        t_app(j,:,i)=[id(j),t_grace];                  
        yhat_app(j,:,i)= [id(j),yihat_grace_raw];
 
 
    end
    
end


%% time onset
id=unique(grace.id);
time_onset=zeros(length(id),4);
group=categorical(length(id),1);

t_begin=nan(length(outcomes),1);
t_end=nan(length(outcomes),1);

for i=1:length(id)
    index_grace=find(grace.id==id(i));
    
    for j=1:length(outcomes)
        tmp_begin=t_app(t_app(:,1,j)==id(i),2,j);
        tmp_end=t_app(t_app(:,1,j)==id(i),3,j);
        if(~isempty(tmp_begin))
            t_begin(j)=tmp_begin;
        end
        if(~isempty(tmp_end))
            t_end(j)=tmp_end;
        end
            
    end
    
    t_min=min(t_begin(~isnan(t_begin)));
    t_max=max(t_end(~isnan(t_end)));   

    time_onset(i,:)=[id(i),t_min,t_max,...
        grace.ConvertTime(index_grace(1))];    
    
    
%     convert=grace.Convert(index_grace(1)); 
%     MCI=grace.MCI(index_grace(1));
%     
%     if(convert)
%          group(i)='pNC';
%     elseif(MCI)
%          group(i)='MCI';
%     else
%          group(i)='sNC';
%     end
    
    convert=grace.Convert(index_grace(1)); 
    
    if(convert==1)
        group(i)='PD-MCI+';
    elseif(convert==0)
        group(i)='PD-MCI-';
    else
        group(i)='NaN';
    end   
    
    t_begin(:)=nan;
    t_end(:)=nan;
    

end



%% Classifier

time_thresh=-2:.2:8;
% acc_thresh_grace=zeros(length(time_thresh),2);
acc_thresh_grace=zeros(length(time_thresh),3);

% MAE_grace=zeros(length(time_thresh),1);


mask_nan= (group=='PD-MCI+') | (group=='PD-MCI-');
mask_Convert= (group=='PD-MCI+')';



for i=1:length(time_thresh)
    %acc_thresh_grace(i,:)=group_classifier_PPMI(time_onset(mask_nan,[2,3]),group(mask_nan),time_thresh(i));
    acc_thresh_grace(i,:)=group_classifier_PPMI_2(time_onset(mask_nan,[2,3]),group(mask_nan),time_thresh(i));
    
    %MAE_grace(i)=mean(abs(time_thresh(i)-time_onset(mask_Convert,2)-time_onset(mask_Convert,4)));


end
% [~,opt_index_grace]=max(sum(acc_thresh_grace,2));
% 
% fprintf('Grace: Time threshold %.1f (range). Accuracy: NC: %.2f, PD: %.2f\n',...
%     time_thresh(opt_index_grace),acc_thresh_grace(opt_index_grace,:)*100);

%%%%%
% fprintf('Grace: Time threshold %.1f (range). Accuracy: PD-CU: %.2f, PD-MCI: %.2f PD: %.2f\n',...
%     time_thresh(opt_index_grace),acc_thresh_grace(opt_index_grace,:)*100);

% [~,opt_index_grace]=min(abs(acc_thresh_grace(:,1)-acc_thresh_grace(:,2)));
[~,opt_index_grace]=min(abs(acc_thresh_grace(:,1)-acc_thresh_grace(:,2))+abs(acc_thresh_grace(:,2)-acc_thresh_grace(:,3)));

fprintf('Grace: Time threshold %.1f (range). Accuracy: PD-CU: %.2f, PD-MCI-BSL: %.2f PD-MCI-END: %.2f\n',...
    time_thresh(opt_index_grace),acc_thresh_grace(opt_index_grace,:)*100);
m_scores=[time_thresh(opt_index_grace),acc_thresh_grace(opt_index_grace,:)*100];
%%%%%
% fprintf('Grace: Time threshold %.1f (range). Accuracy: NC: %.2f, PD: %.2f PD: %.2f\n',...
%     time_thresh(opt_index_grace),acc_thresh_grace(opt_index_grace,:)*100);

% [MAE_grace_min,opt_index_grace]=min(MAE_grace);

% fprintf('Grace: Time threshold %.1f (range). Min MAE: %.2f\n',...
%     time_thresh(opt_index_grace),MAE_grace_min);

%% onset estimation
% time_onset_grace=median(time_onset(group=='PD-MCI+',2)+time_onset(group=='PD-MCI+',4));
% acc_thresh_grace=group_classifier_PPMI(time_onset(:,[2,3]),group,time_onset_grace);

% acc_thresh_grace=group_classifier_PPMI(time_onset(mask_nan,[2,3]),group(mask_nan),time_onset_grace);

% fprintf('Grace: Time threshold %.1f (median). Accuracy: PD-CU: %.2f, PD-MCI: %.2f\n',...
%     time_onset_grace,acc_thresh_grace*100);


%% Correlation
% mask_Convert= (group=='pNC')';
% corr=corrcoef(time_onset_grace-time_onset(mask_Convert,2),time_onset(mask_Convert,4));
% fprintf('Grace: Correlation between onset time and estimation: %.2f\n',corr(1,2));
% m_corr=zeros(2,1);
% m_corr(1)=corr(1,2);
% 
% corr=corrcoef(time_thresh(opt_index_grace)-time_onset(mask_Convert,2),time_onset(mask_Convert,4));
% fprintf('Grace: Correlation between onset time and estimation: %.2f\n',corr(1,2));
% m_corr(2)=corr(1,2);
% 
% 
% if(show_extra)
%     figure(1);clf;
%     subplot(1,2,2);
%     plot(time_onset(mask_Convert,4),time_onset_grace-time_onset(mask_Convert,2),'+')
%     xlabel('Convert time from CU into MCI/dementia (pCU) [years]');
%     ylabel('Estimated convert time, t_{onset}-(t_{i_1}+\delta_i)- [years] ');  
%     title('Grace: estimated convert time (pCU)');
% end

corr=corrcoef(time_thresh(opt_index_grace)-time_onset(mask_Convert,2),time_onset(mask_Convert,4));
fprintf('Grace: Correlation between onset time and estimation (range): %.2f\n',corr(1,2));
%m_corr=zeros(2,1);
m_corr=corr(1,2);

% corr=corrcoef(time_onset_grace-time_onset(mask_Convert,2),time_onset(mask_Convert,4));
% fprintf('Grace: Correlation between onset time and estimation (median): %.2f\n',corr(1,2));
% m_corr(2)=corr(1,2);

% fprintf('Grace: MAE onset time (range): %.2f\n',mean(abs( time_thresh(opt_index_grace)-...
%     time_onset(mask_Convert,2)-time_onset(mask_Convert,4))));
% fprintf('Grace: MAE onset time (median): %.2f\n',mean(abs( time_onset_grace-...
%     time_onset(mask_Convert,2)-time_onset(mask_Convert,4))));
% 
% fprintf('Grace: Rsq (range): %.2f\n',get_Rsq_linear(time_thresh(opt_index_grace)-...
%     time_onset(mask_Convert,2),time_onset(mask_Convert,4)));
% fprintf('Grace: MAE onset time (median): %.2f\n',get_Rsq_linear(time_onset_grace-...
%     time_onset(mask_Convert,2),time_onset(mask_Convert,4)));

% numfigs=length(outcomes);
if(show_extra)
    %figure(numfigs+1);clf;
    figure(1);clf;
    %plot((time_onset_grace-time_onset(mask_Convert,2))*0.3+1.5,time_onset(mask_Convert,4),'+')
    plot(time_onset_grace-time_onset(mask_Convert,2),time_onset(mask_Convert,4),'+')
    xlabel('Convert time estimating');
    ylabel('Convert time from PD-MCI- into PD-MCI+');
    title('Grace');
end

%% Cognitive reservation
if(show_extra)    
    %figure(numfigs+2);clf;
    figure(2);clf;
    histogram(time_onset(group=='PD-MCI-',3)-time_onset_grace);hold on;
    histogram(time_onset(group=='PD-MCI+',3)-time_onset_grace);
    
    xlabel('t_{end}+\gamma_i-t_0 [years]');
    ylabel('Grace: Number of subjects for each \gamma_i');
    legend('PD-MCI-','PD-MCI+');
   
end


%% Longitudinal markers
time_onset_grace=0;
timeline=-10:.1:10;
if(show)
    figure(1);clf;
    num_outcomes=length(outcomes);
%     if(strcmp(cognitiveCriteria{n_cog_criteria},'SITE'))
%         num_outcomes=num_outcomes-1;
%     end
        
    
    for i=1:num_outcomes
        subplot(1,num_outcomes,i);
        %xlabel('years');
        %ylabel('Marker for each subject');
        title(feature_names{i});
        hold on; 
        
        id_tmp=squeeze(t_app(:,1,i));
        id=unique(id_tmp(~isnan(id_tmp)));
 
        t_ghat=[];
        ghat_grace =[];
        for j=1:length(id)
            
            t_grace=t_app(t_app(:,1,i)==id(j),2:3,i)-time_onset_grace;
            yihat_grace_raw=yhat_app(yhat_app(:,1,i)==id(j),2:3,i);

            t_ghat_grace=grace_t_ghat(grace_t_ghat(:,1,i)==id(j),2,i)-time_onset_grace;
            per2feat_ghat_grace=grace_t_ghat(grace_t_ghat(:,1,i)==id(j),3,i);
            t_ghat=[t_ghat;t_ghat_grace];
            ghat_grace=[ghat_grace;per2feat_ghat_grace];

            index=find(grace.id==id(j));
            convert=grace.Convert(index(1));


            if(convert==1)
                plot(t_grace,yihat_grace_raw,'r');
                %plot(t_ghat_grace,per2feat_ghat_grace,'ro');
            elseif(convert==0)
                plot(t_grace,yihat_grace_raw,'g');
                %plot(t_ghat_grace,per2feat_ghat_grace,'g+');
            else
            end           
 
        end
        %plot(t_ghat,ghat_grace,'k+');
        [t_ghat,idx]=sort(t_ghat);
        ghat_grace=ghat_grace(idx);
        [t_ghat,idx,~]=unique(t_ghat);
        ghat_grace=ghat_grace(idx);
        plot(timeline,interp1(t_ghat,ghat_grace,timeline,'pchip'),...
             '-','Color','k','LineWidth',2);
        xlim([-17 17]);
        if(original_scale==false)
            ylim([0 1]);
        end
       
       
       
        grid on;

    end
end






    
%% exit
%rmpath('./aux_grace');

end
