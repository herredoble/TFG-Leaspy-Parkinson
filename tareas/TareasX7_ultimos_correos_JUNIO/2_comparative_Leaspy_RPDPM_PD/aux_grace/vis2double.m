function years = vis2double(vis)
years=zeros(length(vis),1);
for i=1:length(vis)
    years(i)= viscode2double(vis(i));
end

end



function visit_year = viscode2double(vis)
if(vis=='BL')
    visit_year=0;
elseif(vis=='V04')
    visit_year=1;
elseif(vis=='V06')
    visit_year=2;
elseif(vis=='V08')
    visit_year=3;
elseif(vis=='V10')
    visit_year=4;
elseif(vis=='V12')
    visit_year=5;
elseif(vis=='ST') %3564
    visit_year=1;
else
    visit_year=NaN;
end
 
end

