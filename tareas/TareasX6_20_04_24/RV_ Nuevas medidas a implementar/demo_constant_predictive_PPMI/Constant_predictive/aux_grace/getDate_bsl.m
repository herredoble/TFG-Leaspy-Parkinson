function date_bsl=getDate_bsl(PATNO,date,YEAR)
date_bsl=date;
RID=unique(PATNO);
for i=1:length(RID)
    date_start=date(PATNO==RID(i) & YEAR==0);
    date_bsl(PATNO==RID(i))=date_start;
end

end