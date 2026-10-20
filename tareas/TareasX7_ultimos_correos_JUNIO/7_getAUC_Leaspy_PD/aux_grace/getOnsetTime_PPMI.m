function onset=getOnsetTime_PPMI(ID,date,convert,convertTime,YEAR)
onset=nan(length(ID),1);
%% Convert
RID=unique(ID(convert==1));
for i=1:length(RID)
    index=find(ID==RID(i));
    convertTime_year=convertTime(index(1));
    years=YEAR(index);
    idx=find(years==convertTime_year,1);
    time_zero=date(index(idx));
    try
        Delta_days=days(date(index(idx))-date(index(idx-1)))/2;
        onset(index)=(days(date(index)-time_zero)+Delta_days)/30.4167/12;
    catch   
        onset(index)=(days(date(index)-date(index(1))))/30.4167/12;
    end
        
end



%% sCU (clinical) 419s 2142v
RID=unique(ID(convert==0));
for i=1:length(RID)
    index=find(ID==RID(i));
    time_zero=date(index(end));
    onset(index)=days(date(index)-time_zero)/30.4167/12;    
    
end

end
