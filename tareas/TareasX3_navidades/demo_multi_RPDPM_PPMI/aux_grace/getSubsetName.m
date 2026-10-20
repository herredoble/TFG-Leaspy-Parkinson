function name = getSubsetName(features_subset_table, i, head_name)
    name = head_name;
    for j = 1 : features_subset_table{i,'Size'}
        name = strcat(name,'_',string(features_subset_table{i,'Features_Subset'}(j)));
    end
end