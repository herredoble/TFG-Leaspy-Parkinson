function ADNI=newADNIMERGE(n_markers)
type_markers={'MRI+NM','MRI+NM+CSF'};


% %addpath('./SelectedVisit_aux');
% %% Change VISCODE,COLPROT,ORIGPR,FSVERSION
% %% Empty FLDSTRENG, FLDSTRENG_bsl
% % %load('./data/ADNIMERGE_201019.mat','ADNIMERGE');
% % load('./data/ADNIMERGE_201118.mat','ADNIMERGE');
% % 
% % ADNI=ADNIMERGE;
% % clear ADNIMERGE
% 
% %% ADNIMERGE R
% ADNI=readtable('./data/ADNIMERGE_R_201123.xlsx');
% ADNI.COLPROT=categorical(ADNI.COLPROT);
% ADNI.ORIGPROT=categorical(ADNI.ORIGPROT);
% ADNI.PTID=categorical(ADNI.PTID);
% ADNI.SITE=categorical(ADNI.SITE);
% ADNI.VISCODE=categorical(ADNI.VISCODE);
% ADNI.DX_bl=categorical(ADNI.DX_bl);
% ADNI.PTGENDER=categorical(ADNI.PTGENDER);
% ADNI.PTETHCAT=categorical(ADNI.PTETHCAT);
% ADNI.PTRACCAT=categorical(ADNI.PTRACCAT);
% ADNI.PTMARRY=categorical(ADNI.PTMARRY);
% 
% ADNI.FLDSTRENG=categorical(ADNI.FLDSTRENG);
% ADNI.FSVERSION=categorical(ADNI.FSVERSION);
% ADNI.DX=categorical(ADNI.DX);
% ADNI.FLDSTRENG_bl=categorical(ADNI.FLDSTRENG_bl);
% ADNI.FSVERSION_bl=categorical(ADNI.FSVERSION_bl);
% 
% ADNI.Month=str2double(ADNI.Month);
% %% Problem with RID 6014 and 0 visit.
% ADNI(ADNI.RID==6014,:)=[];
% 
% 
% %% info
% clc;
% 
% subjects=unique(ADNI.RID);
% fprintf('ADNIMERGE new version: %d subjects with %d visits.\n',...
%     length(subjects),size(ADNI,1));
% 
% %% Replace FLDSTRENG with data from ADNIMERGE_14036v_2175s
% % load('./data/ADNIMERGE_14036v_2175s.mat','ADNIMERGE_14036v_2175s');
% % subjects_old=unique(ADNIMERGE_14036v_2175s.RID);
% % fprintf('ADNIMERGE old version: %d subjects with %d visits.\n',...
% %     length(subjects_old),size(ADNIMERGE_14036v_2175s,1));
% % fprintf('Change FLDSTRENG values from old to new version.\n');
% % pause;
% % 
% % addpath('./SelectedVisit_aux');
% % index=changeVariableADNIMERGE(ADNI,ADNIMERGE_14036v_2175s);
% % FLDSTRENG=ADNIMERGE_14036v_2175s.FLDSTRENG(index(index>0));
% % mask=(FLDSTRENG=='1.5 Tesla MRI') | (FLDSTRENG=='3 Tesla MRI');
% % FLDSTRENG(~mask)='NaN';
% % ADNI.FLDSTRENG=categorical(ADNI.FLDSTRENG);
% % ADNI.FLDSTRENG_bl=categorical(ADNI.FLDSTRENG_bl);
% % 
% % ADNI.FLDSTRENG(index>0)=FLDSTRENG;
% % ADNI.FLDSTRENG(index==0)='NaN';
% % 
% % fprintf('New visits without FLDSTRENG data: %d, Total new visit: %d.\n',...
% %     sum(index==0),size(ADNI,1)-size(ADNIMERGE_14036v_2175s,1));
% 
% %% FLDSTRENG_bl
% % for i=1:length(subjects)
% %     index=find(ADNI.RID==subjects(i));
% %     FDL_bsl=ADNI.FLDSTRENG(index(ADNI.M(index)==0));
% %     try
% %         ADNI.FLDSTRENG_bl(index)=FDL_bsl;
% %     catch
% %         ADNI.FLDSTRENG_bl(index)='NaN';
% %     end
% % end
% 
% %% Checking
% % fprintf('Finished task. Now checking')
% % pause;
% % clc;
% % index=changeVariableADNIMERGE(ADNIMERGE_14036v_2175s,ADNI);
% % if(sum(index==0))
% %     fprintf('There are %d visits from old version without reference in new version.\n',sum(index==0));
% % end
% % 
% % FLD=[ADNIMERGE_14036v_2175s.FLDSTRENG(index>0),ADNI.FLDSTRENG(index(index>0)),...
% %      ADNIMERGE_14036v_2175s.FLDSTRENG_bl(index>0),ADNI.FLDSTRENG_bl(index(index>0))];
% % mask=(FLD=='1.5 Tesla MRI') | (FLD=='3 Tesla MRI');
% % FLD(~mask)='NaN';
% % if(sum(FLD(:,1)~=FLD(:,2))>0)
% %     fprintf('Discrepances between old and new version from %d visits in FLDSTRENG data',...
% %         sum(FLD(:,2)~=FLD(:,3)));
% %         
% % end
% % 
% % if(sum(FLD(:,3)~=FLD(:,4))>0)
% %     fprintf('Discrepances between old and new version from %d visits in FLDSTRENG baseline data',...
% %         sum(FLD(:,2)~=FLD(:,3)));
% %         
% % end
% %% CSF
% ABETA=str2double(ADNI.ABETA);
% char_ABETA=char(ADNI.ABETA);
% mask_Gr=strcmp(string(char_ABETA(:,1)),'>');
% ABETA(mask_Gr)=1700;
% mask_ls=strcmp(string(char_ABETA(:,1)),'<');
% ABETA(mask_ls)=200;
% ADNI.ABETA=ABETA;
% 
% TAU=str2double(ADNI.TAU);
% char_TAU=char(ADNI.TAU);
% mask_Gr=strcmp(string(char_TAU(:,1)),'>');
% TAU(mask_Gr)=1300;
% mask_ls=strcmp(string(char_TAU(:,1)),'<');
% TAU(mask_ls)=80;
% ADNI.TAU=TAU;
% 
% PTAU=str2double(ADNI.PTAU);
% char_PTAU=char(ADNI.PTAU);
% mask_Gr=strcmp(string(char_PTAU(:,1)),'>');
% PTAU(mask_Gr)=120;
% mask_ls=strcmp(string(char_PTAU(:,1)),'<');
% PTAU(mask_ls)=8;
% ADNI.PTAU=PTAU;
% %% Baseline summary
% % fprintf('verification completed. Baseline summary.\n');
% % pause;
% % clc;
% % ADNI_bsl=ADNI(ADNI.M==0,:);
% % features=[ADNI_bsl.FAQ,ADNI_bsl.RAVLT_immediate,abs(ADNI_bsl.RAVLT_perc_forgetting),...
% %     ADNI_bsl.ADASQ4,ADNI_bsl.Hippocampus./ADNI_bsl.ICV*1e3];
% % feature_names={'FAQ','RAVLT_immediate','RAVLT_perc_forgetting','ADASQ4','NHV'};
% % summary_bsl(features,ADNI_bsl.DX,feature_names);
% 
% %% Error in markers
% % fprintf('Number of visits where RAVLTpercForgetting is more than 100: %d\n',...
% %     sum(abs(ADNI.RAVLT_perc_forgetting)>100));
% % ADNI.RAVLT_perc_forgetting(abs(ADNI.RAVLT_perc_forgetting)>100)=nan;
% % ADNI.RAVLT_perc_forgetting_bl(abs(ADNI.RAVLT_perc_forgetting_bl)>100)=nan;
% % 
% % ADNI_bsl=ADNI(ADNI.M==0,:);
% % features=abs(ADNI_bsl.RAVLT_perc_forgetting);
% % feature_names={'RAVLT_perc_forgetting'};
% % summary_bsl(features,ADNI_bsl.DX,feature_names);


%% Filters
if(strcmp(type_markers{n_markers},'MRI+NM'))
    %ADNI = mask_MRI_NM(ADNI);
    load('./data/ADNIMERGE_MRI_NM_201202','ADNI');
elseif(strcmp(type_markers{n_markers},'MRI+NM+CSF'))
    %ADNI = mask_MRI_NM_CSF(ADNI);
    load('./data/ADNIMERGE_MRI_NM_CSF_201202','ADNI');
else
end

%% delete
%rmpath('./SelectedVisit_aux');


end