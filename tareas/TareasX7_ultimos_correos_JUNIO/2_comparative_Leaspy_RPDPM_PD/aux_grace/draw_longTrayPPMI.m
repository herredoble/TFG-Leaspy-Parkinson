function draw_longTrayPPMI(RID,time,feat,group,level_linear)
uRID=unique(RID);
for i=1:length(uRID)
    if(sum(RID==uRID(i))>=2)
        grp=unique(group(RID==uRID(i)));
        drawSplines(feat(RID==uRID(i)),time(RID==uRID(i)),level_linear,grp);
    end
end

end

function drawSplines(samples,time,lineal,group)
if(strcmp(group,'PD'))
    index=2;
elseif(strcmp(group,'Control'))
    index=1;
else
    index=3;
end
colors={'g','b','r','k'};
sp=csaps(time,samples,lineal);
fnplt(sp,colors{index});

end

