function timeline=get_timeline_tray(ages,years_before,years_after,resolution_year)

J2=max_samples_trayectories(ages,years_before,years_after,resolution_year);
timeline= NaN(size(ages,1),J2);
for i=1:size(ages,1)
    t_max=max(ages(i,:));
    t_min=min(ages(i,:));
    samples=t_min-years_before:resolution_year:t_max+years_after;
    timeline(i,1:length(samples))=samples;
end

end

function n_tray_samples=max_samples_trayectories(ages,years_before,years_after,resolution_year)

n_tray_samples=0;
for i=1:size(ages,1)
    t_max=max(ages(i,:));
    t_min=min(ages(i,:));
    samples=t_min-years_before:resolution_year:t_max+years_after;
    if(length(samples)>n_tray_samples)
        n_tray_samples=length(samples);
    end
end
end