function [feat2per,percent,point_feat]=percentil_interp(features,group)
%function [feat2per,percent,point_feat]=percentil_interp(features,group,percent,point_feat)

numFeatures=size(features,2);
feat2per=zeros(size(features));
if(nargin>2)
    for i=1:numFeatures
        feat2per(:,i)=interp1(point_feat{i},percent{i},features(:,i));     
    end
else
    percent=cell(numFeatures,1);
    point_feat=cell(numFeatures,1);
    for i=1:numFeatures
        mask_nan=~isnan(features(:,i));
        [prob,feat]=ecdf(balanceSamples(features(mask_nan,i),group(mask_nan)));
%         point_feat{i}=feat(2:end);
%         percent{i}=prob(2:end);
        feat(2)=feat(1)+(feat(3)-feat(2))/2;
        point_feat{i}=feat;
        percent{i}=prob;
        
        feat2per(:,i)=interp1(point_feat{i},percent{i},features(:,i));     
    end
    
end


end

function features=balanceSamples(features,group)
gr_name=unique(group);
num_samples=zeros(length(gr_name),1);
for i=1:length(gr_name)
    num_samples(i)=sum(group==gr_name(i));
end
max_samples=max(num_samples);
feat_extra=[];
for i=1:length(gr_name)
    if((max_samples-num_samples(i))>0)
       feat=features(group==gr_name(i));
       index=randi(num_samples(i),max_samples-num_samples(i),1);
       feat_extra=[feat_extra;feat(index)];
    end
end

features=[features;feat_extra];
        

end
