function [LMEmodels,MAE]=getLME_univariate_corr(Tresults,features)

num_features=size(features,2);
LMEmodels=cell(num_features,1);


years=Tresults.Month/12;

mask_bsl=years==0;
MAE=nan(2*sum(mask_bsl),num_features);

for i=1:num_features
    %% LME
    mask_noNaN=getMasknoNaNnoBsl(features(:,i),mask_bsl,Tresults.RID);
    feat_long=features(mask_noNaN,i);
    years_long=years(mask_noNaN);
    numSubjs=numel(unique(Tresults.RID(mask_noNaN)));
       
    ni=get_ni(years_long,numSubjs);
    LMEmodels{i}  = getLME_Model(feat_long,years_long,ni);
    %% MAE
    for j=1:length(ni)
        if(j>1)
            index_begin=index_end+1;
        else
            index_begin=1;        
        end
        index_end=index_begin+ni(j)-1; 
        years_ID=years_long(index_begin:index_end);
        [tmin,idx_min]=min(years_ID);
        [tmax,idx_max]=max(years_ID);

        feat_samples=feat_long(index_begin:index_end);

        y_bsl=LMEmodels{i}.Bhat(1)+(LMEmodels{i}.bihat(1,j));
        y_end=y_bsl+ ((tmax-tmin)*(LMEmodels{i}.Bhat(2)+LMEmodels{i}.bihat(2,j)));
        MAE(j*2-1:j*2,i)=[abs(feat_samples(idx_min)-y_bsl);....
            abs(feat_samples(idx_max)-y_end)];

    end


end

end

function model  = getLME_Model(marker,time_visit,ni)

numImgs=length(marker);
X_Hipp=[ones(numImgs,1),time_visit];
model = lme_fit_FS(X_Hipp,[1 2],marker,ni);


end


function mask_noNaN=getMasknoNaNnoBsl(feature,mask_bsl,RID)

mask_noNaN=~isnan(feature);
index_subjs_withoutBsl = find(mask_bsl & mask_noNaN==0);
for i=1:length(index_subjs_withoutBsl)
    mask_noNaN(RID==RID(index_subjs_withoutBsl(i)))=false;
end


end


