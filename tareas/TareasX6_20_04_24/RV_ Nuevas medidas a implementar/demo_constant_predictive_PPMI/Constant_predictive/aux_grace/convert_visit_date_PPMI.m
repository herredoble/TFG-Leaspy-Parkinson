function date=convert_visit_date_PPMI(visit_date)

n_visit=length(visit_date);
visit_date=char(visit_date);
%month
month_str=string(visit_date(:,1:3));
month=zeros(n_visit,1);
month(month_str=="JAN")=1;
month(month_str=="FEB")=2;
month(month_str=="MAR")=3;
month(month_str=="APR")=4;
month(month_str=="MAY")=5;
month(month_str=="JUN")=6;
month(month_str=="JUL")=7;
month(month_str=="AUG")=8;
month(month_str=="SEP")=9;
month(month_str=="OCT")=10;
month(month_str=="NOV")=11;
month(month_str=="DEC")=12;
%year
year=zeros(n_visit,1);
for i=1:n_visit    
    year(i)=str2double(visit_date(i,4:end));   
end

%date
date = arrayfun(@(m, y) datetime(y, m, 1), month, year);
end