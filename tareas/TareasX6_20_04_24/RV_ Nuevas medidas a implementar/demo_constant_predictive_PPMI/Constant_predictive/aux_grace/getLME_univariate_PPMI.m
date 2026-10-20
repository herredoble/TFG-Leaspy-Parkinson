function LMEmodels=getLME_univariate_PPMI(Tresults,features,feature_names,...
    bMRImarker,SxEdGr,bTwoGroup)

num_features=size(features,2);
LMEmodels=cell(num_features,1);


% years=Tresults.Month/12;
years=Tresults.years;
mask_bsl=years==0;

%Age=Tresults.Age-years;
Age=Tresults.age;

% meanAge=mean(Age(mask_bsl));
% meanEdu=mean(Tresults.EDUCYRS(mask_bsl));
% meanAge=61.57;
% meanEdu=15.84;
meanAge=60.5;
meanEdu=15.59;

 

Age=Age-meanAge;
Education=Tresults.EDUCYRS-meanEdu;
% Sex=double(strcmp(Tresults.sex,'M'));
Sex=Tresults.gen;



for i=1:num_features
    fprintf('Model %s\n',feature_names{i});
    mask_noNaN=getMasknoNaNnoBsl(features(:,i),mask_bsl,Tresults.ID);
    feat_long=features(mask_noNaN,i);
    years_long=years(mask_noNaN);
    Age_long=Age(mask_noNaN);
    fld_long=[];
    numSubjs=numel(unique(Tresults.ID(mask_noNaN)));
    
    if(SxEdGr)
        Sex_long=Sex(mask_noNaN);
        Edu_long=Education(mask_noNaN);
    else
        Sex_long=[];
        Edu_long=[];
    end
    
    if(bTwoGroup)
        Grp_long=Convert(mask_noNaN);
    else
        Grp_long=[];
    end
    
    ni=get_ni(years_long,numSubjs);
    if(bMRImarker(i))
        LMEmodels{i}  = getLME_Model_field(feat_long,fld_long,years_long,...
            Age_long,Grp_long,Sex_long,Edu_long,ni);
    else
        LMEmodels{i}  = getLME_Model(feat_long,years_long,Age_long,...
            Grp_long,Sex_long,Edu_long,ni);
    end

end

end

function model  = getLME_Model(marker,time_visit,blAge,convert,Sex,Education,ni)

numImgs=length(marker);
if(isempty(convert) && isempty(Sex))
    X_Hipp=[ones(numImgs,1),time_visit,blAge];
elseif(isempty(Sex))
    X_Hipp=[ones(numImgs,1),time_visit,blAge,convert,convert.*time_visit];
elseif(isempty(convert))
    X_Hipp=[ones(numImgs,1),time_visit,blAge,Sex,Education];
else
    X_Hipp=[ones(numImgs,1),time_visit,blAge,convert,convert.*time_visit,Sex,Education];

end
model = lme_fit_FS(X_Hipp,[1 2],marker,ni);


end


function model  = getLME_Model_field(marker,field,time_visit,blAge,convert,Sex,Education,ni)

numImgs=length(marker);
Teslas=unique(field);

if(isempty(convert) && isempty(Sex))
    X_Hipp=[ones(numImgs,1),time_visit,blAge,double(field==Teslas(end))];
elseif(isempty(Sex))
    X_Hipp=[ones(numImgs,1),time_visit,blAge,double(field==Teslas(end)),...
        convert,convert.*time_visit];
elseif(isempty(convert))
    X_Hipp=[ones(numImgs,1),time_visit,blAge,double(field==Teslas(end)),...
        Sex,Education];
else
    X_Hipp=[ones(numImgs,1),time_visit,blAge,double(field==Teslas(end)),...
        convert,convert.*time_visit,Sex,Education];    
end


model = lme_fit_FS(X_Hipp,[1 2],marker,ni);


end

function mask_noNaN=getMasknoNaNnoBsl(feature,mask_bsl,RID)

mask_noNaN=~isnan(feature);
index_subjs_withoutBsl = find(mask_bsl & mask_noNaN==0);
for i=1:length(index_subjs_withoutBsl)
    mask_noNaN(RID==RID(index_subjs_withoutBsl(i)))=false;
end


end


