function features_corr=correctionCovariablesLME_PPMI(LMEmodels,Tresults,features,ni)

numSubjects=length(ni);
numFeatures=numel(LMEmodels);

% years=Tresults.years;
% age=Tresults.Age-years;
age=Tresults.age;

% sex=double(strcmp(Tresults.sex,'M'));
sex=Tresults.gen;

features_corr=NaN(size(features));

% meanAge=61.57;
% meanEdu=15.84;
meanAge=60.5;
meanEdu=15.59;

for i=1:numSubjects
    if(i>1)
        index_begin=index_end+1;
    else
        index_begin=1;        
    end
    index_end=index_begin+ni(i)-1;
    %Covariables
    Age=age(index_begin)-meanAge;
    Sex=sex(index_begin);
    Edu=Tresults.EDUCYRS(index_begin)-meanEdu;
    
    %Var
    feat_samples=features(index_begin:index_end,:);
    feat_estimated=nan(ni(i),numFeatures);
    for j=1:numFeatures
        mask_nan=isnan(feat_samples(:,j))==0;
        gi=getCorrMarkerTrayectorySubject(LMEmodels{j},Age,[],...
            Sex,Edu,feat_samples(mask_nan,j));      
        feat_estimated(mask_nan,j)=gi;
    end
    
   features_corr(index_begin:index_end,:)=feat_estimated;
    
    
end


end




function gi=getCorrMarkerTrayectorySubject(model,Age,fld,sex,edu,yi)
% if(length(model.Bhat)==4)%X_Hipp=[ones,time_visit,blAge,field==Teslas(end)];
%     %bsl=model.Bhat(4)*Age);
%     bsl=model.Bhat(3)*Age;  
%     gi= yi-repmat(bsl,length(yi),1)-(model.Bhat(4)*fld);
    
if(length(model.Bhat)==3)%X_Hipp=[ones,time_visit,blAge];
    bsl=model.Bhat(3)*Age;
    gi= yi-repmat(bsl,length(yi),1);

elseif(length(model.Bhat)==5)%X_Hipp=[ones(numImgs,1),time_visit,blAge,Sex,Education];
    bsl=(model.Bhat(3)*Age)+(model.Bhat(4)*sex)+(model.Bhat(5)*edu);
    gi= yi-repmat(bsl,length(Age),1);
    
% elseif(length(model.Bhat)==7)%X_Hipp=[ones,time_visit,blAge,convert,convert*time,sex,edu];
%     bsl=(model.Bhat(3)*Age)+(model.Bhat(6)*sex)+(model.Bhat(7)*edu);
%     gi= yi-repmat(bsl,length(Age),1);
%         
% elseif(length(model.Bhat)==8)%X_Hipp=[ones,time_visit,blAge,convert,convert*time,sex,edu];
%     bsl=(model.Bhat(3)*Age)+(model.Bhat(7)*sex)+(model.Bhat(8)*edu);
%     gi= yi-repmat(bsl,length(Age),1)-(model.Bhat(4)*fld); 
    
else
    gi=[];
end
end