function [inter,slope]=corrOutcomeLME(LMEmodels)
numOut=numel(LMEmodels);
numSub=length(LMEmodels{1}.ni);
inter=zeros(numSub,numOut);
slope=zeros(numSub,numOut);
for i=1:numOut
    if(length(LMEmodels{i}.Bhat)==2)
        inter(:,i)=LMEmodels{i}.Bhat(1)+LMEmodels{i}.bihat(1,:)';
        slope(:,i)=LMEmodels{i}.Bhat(2)+LMEmodels{i}.bihat(2,:)';
    elseif(length(LMEmodels{i}.Bhat)==3) % Age(bsl)
        inter(:,i)=LMEmodels{i}.bihat(1,:)';
        mask_bsl=LMEmodels{i}.X(:,2)==0;
        inter(:,i)=inter(:,i)+ (LMEmodels{i}.X(mask_bsl,3)*LMEmodels{i}.Bhat(3));
        inter(:,i)=inter(:,i)+  LMEmodels{i}.Bhat(1);
        slope(:,i)=LMEmodels{i}.bihat(2,:)';
        slope(:,i)=slope(:,i)+LMEmodels{i}.Bhat(2);
    elseif(length(LMEmodels{i}.Bhat)==4) % Age(bsl)+FDL(bsl)
        inter(:,i)=LMEmodels{i}.bihat(1,:)';
        mask_bsl=LMEmodels{i}.X(:,2)==0;
        inter(:,i)=inter(:,i)+ (LMEmodels{i}.X(mask_bsl,3)*LMEmodels{i}.Bhat(3))+...
            (LMEmodels{i}.X(mask_bsl,4)*LMEmodels{i}.Bhat(4));
        inter(:,i)=inter(:,i)+  LMEmodels{i}.Bhat(1);        
        slope(:,i)=LMEmodels{i}.bihat(2,:)'; 
        slope(:,i)=slope(:,i)+LMEmodels{i}.Bhat(2);        
    end
end

%corr_inter=corrcov(cov(inter));
%corr_slope=corrcov(cov(slope));

end
