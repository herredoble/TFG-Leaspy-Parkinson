function visit=visit2double(VISCODE)

visit=-ones(length(VISCODE),1);
for i=1:length(VISCODE)
    visit(i)=visit2num(VISCODE(i));
end

end


function visit=visit2num(VISCODE)
vis=char(VISCODE);
if(vis(1)=='b')
    visit=0;
else
    visit=str2double(vis(2:end));
end

end