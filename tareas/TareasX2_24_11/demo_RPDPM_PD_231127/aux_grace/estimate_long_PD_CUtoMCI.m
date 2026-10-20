function [raw_feat,feature_names,PPMI,long_term]=estimate_long_PD_CUtoMCI(original_scale,t_app)

addpath('./aux_grace');


%% data
n_dataset=1;
n_cog_criteria=3;


num_outcomes=size(t_app,3);

n_typeFeat=5;
[PPMI,~,~,~,~] = getDataPPMI_MCI_v1(n_typeFeat,n_dataset);
% Cognitive criteria
cognitiveCriteria={'MCI_test','UPDRS1','SITE'};
varCogCriteria=cellstr(cognitiveCriteria{n_cog_criteria});
if(strcmp(varCogCriteria,'SITE'))
    PPMI.Convert=double(PPMI.cogstate_end>1);
    PPMI.convertTime=PPMI.convertTimeCogstate;
elseif(strcmp(varCogCriteria,'MCI_test'))
    PPMI.Convert=double(PPMI.MCI_testscores_end);
    PPMI.convertTime=PPMI.convertTimeMCI_test;
elseif(strcmp(varCogCriteria,'UPDRS1'))
    PPMI.Convert=double(PPMI.np1cog_end>1);
    PPMI.convertTime=PPMI.convertTimeNp1cog;
else
end

PPMI.Month=PPMI.YEAR*12;



%% feature add
RID=unique(PPMI.PATNO);

% feature=[PPMI.moca,PPMI.abeta,PPMI.ptau,PPMI.ptau_ab];
% feature_names=["MOCA","A\beta","pTAU","pTAU/A\beta"];
% change_sign_feat=[true;true;false;false];

% feature=[PPMI.SDMTOTAL,PPMI.HVLTRDLY,PPMI.sft,PPMI.lns,...
%     PPMI.moca,PPMI.mean_caudate,PPMI.ptau_ab,PPMI.abeta];
% feature_names=["SDMT","HVLTRDLY","SFT","LNS","MOCA","MEANCAU","pTAU/A\beta","A\beta"];
% change_sign_feat=true(length(feature_names),1);
% change_sign_feat(end-1)=false;
% show =true;



% feature=[PPMI.SDMTOTAL,PPMI.lns,PPMI.abeta];
% feature_names=["SDMT","LNS","A\beta"];
% change_sign_feat=true(length(feature_names),1);
% show=false;

% feature=[PPMI.moca,PPMI.abeta,PPMI.SDMTOTAL,PPMI.HVLTRDLY,PPMI.sft,...
%     PPMI.mean_caudate,PPMI.ptau];
% feature_names=["MOCA","A\beta","SDMT","HVLTRDLY","SFT","MEANCAUDATE","pTAU"];
% change_sign_feat=true(length(feature_names),1);
% change_sign_feat(end)=false;

feature=[PPMI.moca,PPMI.abeta,PPMI.SDMTOTAL,PPMI.mean_caudate];
feature_names=["MOCA","A\beta","SDMT","MEANCAUDATE"];
change_sign_feat=true(length(feature_names),1);

% feature=[PPMI.moca,PPMI.lns,PPMI.ptau];
% feature_names=["MOCA","LNS","pTAU"];
% change_sign_feat=true(length(feature_names),1);
% change_sign_feat(end)=false;

%% Transforming
if(~original_scale)
    for i=1:length(change_sign_feat)
        if(change_sign_feat(i)==true)
            feature(:,i)=-feature(:,i);
        end
    end    
    [feature,~,~]=percentil_interp(feature,PPMI.Convert);
end


%% Features
% num_fig=1;
% figure(num_fig);clf;

t_begin=nan(num_outcomes,1);
timeline=-15:.1:15;
long_term=zeros(length(timeline),size(feature,2)+1);
long_term(:,1)=timeline';
raw_feat=[];

for k=1:size(feature,2)
    samples_=nan(size(feature,1),4);
    idx=1;
    
for i=1:length(RID)
    mask_subj=PPMI.PATNO==RID(i);
    
    feat_subj=feature(mask_subj,k);
    mask_nan=~isnan(feat_subj);
    feat_subj=feat_subj(mask_nan);
    
    if(~isempty(feat_subj))
    
        PPMI_subj=PPMI(PPMI.PATNO==RID(i),:);
        if(issorted(PPMI_subj.YEAR))
            t_subj=PPMI_subj.YEAR(mask_nan);
            
            t_begin=nan(num_outcomes,1);
            for j=1:num_outcomes
                tmp_begin=t_app(t_app(:,1,j)==RID(i),2,j);
                if(~isempty(tmp_begin))
                    t_begin(j)=tmp_begin;
                end
            end
            t_offset=min(t_begin(~isnan(t_begin)));
            %t_begin(:)=nan;

            num_samples=length(feat_subj);
            
            %%%%%%%%%%%%%%%%%%%%
%             if(num_samples>2)
%                 %feat_subj_int=interp1(t_subj+t_offset,feat_subj,t_subj+t_offset);
%                 [pp,SS] = polyfit(t_subj+t_offset,feat_subj,1);
%                 [feat_subj_int,~] = polyval(pp,t_subj+t_offset,SS);
%             else
                feat_subj_int=feat_subj;
%             end
            %%%%%%%%%%%%%%%%%%%%
            
            samples_(idx:idx+num_samples-1,:)=[repmat(RID(i),num_samples,1),...
                t_subj+t_offset,feat_subj_int,repmat(PPMI_subj.Convert(1),num_samples,1)];
            idx=idx+num_samples;

%             if(show)
%                 t_feat=[t_subj(1),t_subj(end)]+t_offset;
%                 if(length(feat_subj)>2)
%                     v_feat=linear_reg(t_subj,feat_subj);
%                     v_feat=[v_feat(1),v_feat(end)];
% 
%                 else
%                     v_feat=[feat_subj(1),feat_subj(end)];
%                 end
% 
%                 if(strcmp(varCogCriteria,'SITE'))
%                     if(PPMI_subj.cogstate_end(1)==1)
%                         plot(t_feat,v_feat,'g','LineWidth',1);
%                     elseif(PPMI_subj.cogstate_end(1)==2) 
%                         plot(t_feat,v_feat,'b','LineWidth',1);
%                     elseif(PPMI_subj.cogstate_end(1)==3)
%                         plot(t_feat,v_feat,'r','LineWidth',1);
%                     else
%                     end
%                 end
%             end

        else
            fprintf(2,'RID %d is not sorted years\n',RID(i));
        end
    else
        fprintf(2,'RID %d is empty for feature %d\n',RID(i),k);
    end

end
samples_=samples_(~isnan(samples_(:,1)),:);

percentiles=50;
q=percentiles/100;
[p,~]=quantreg(samples_(:,2),samples_(:,3),q,2);
long_term(:,k+1)=polyval(p,timeline');
raw_feat=[raw_feat;samples_(:,1),samples_(:,2),samples_(:,3),repmat(k,size(samples_,1),1)];

end


%% exit
% rmpath('./aux_grace');
end


