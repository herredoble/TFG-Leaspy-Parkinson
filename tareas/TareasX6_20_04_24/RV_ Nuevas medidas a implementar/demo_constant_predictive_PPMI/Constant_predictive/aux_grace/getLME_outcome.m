function model=getLME_outcome(ADNI)



ID=unique(ADNI.id);
ni=zeros(length(ID),1);
t_study=zeros(size(ADNI,1),1);
for i=1:length(ID)
    index=find(ADNI.id==ID(i));
    ni(i)=length(index);
    tmin=min(ADNI.t(index));
    t_study(index)=ADNI.t(index)-tmin;   
end
model  = getLME_Model(ADNI.Y,t_study,ni);


end







function model  = getLME_Model(marker,time_visit,ni)

numImgs=length(marker);
X_Hipp=[ones(numImgs,1),time_visit];
model = lme_fit_FS(X_Hipp,[1 2],marker,ni);


end