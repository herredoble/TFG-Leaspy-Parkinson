function grace=cleanID_sinVisit(grace)
outcomes=unique(grace.Outcome);
IDs=unique(grace.id);
m_ID=nan(length(IDs),length(outcomes));
for i=1:length(outcomes)
    ID=unique(grace.id(grace.Outcome==outcomes(i)));
    m_ID(1:length(ID),i)=ID;
end
[C,~,ic] = unique(m_ID);
hist_ID = [C, accumarray(ic,1)];
del_ID=hist_ID(hist_ID(:,2)~=length(outcomes),1);
del_ID=del_ID(~isnan(del_ID));
index=[];
for i=1:length(del_ID)
    index=[index;find(grace.id==del_ID(i))];
end

grace(index,:)=[];

end