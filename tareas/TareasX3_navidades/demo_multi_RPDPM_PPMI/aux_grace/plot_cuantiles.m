function plot_cuantiles
load('./data/poblacion_pNC','final_pNC_onset');
data=table;
data.x=final_pNC_onset.onset;
% data.y=final_pNC_onset.ADAS13;
% data.y=final_pNC_onset.MMSE;
data.y=final_pNC_onset.FAQ;

% NHV=final_pNC_onset.Hippocampus./final_pNC_onset.ICV*1e3;
% NHV=final_pNC_onset.Entorhinal./final_pNC_onset.ICV*1e3;
% NHV=final_pNC_onset.Ventricles./final_pNC_onset.ICV*1e3;
% NHV=final_pNC_onset.MidTemp./final_pNC_onset.ICV*1e3;
% NHV=final_pNC_onset.Fusiform./final_pNC_onset.ICV*1e3;

% data.y=NHV;

data=sortrows(data);
x=data.x(~isnan(data.y));
y=data.y(~isnan(data.y));
% percentiles=[5 25 50 75 95];
percentiles=[50-34.1 50 50+34.1];

q=percentiles/100;
figure(1);clf;
hold on
plot(x,y,'+');
for i=1:length(q)
    [p,stats]=quantreg(x,y,q(i),2);
    plot(x,polyval(p,x))
%     if(i==3)
%         plot(x,y,x,polyval(p,x),x,stats.yfitci,'k:')
%     else 
%         plot(x,polyval(p,x))
%     end
end
% legend('data','5th percentile','25th percentile','50th percentile','75th percentile','95th percentile')
end
