function [fd_summ] = get_dataset_summary(metrics)
    fields = fieldnames(metrics);

    for i = 1 : length(fields) 
        if fields{i} ~= "File"
            tempVect = [metrics(:).(fields{i})];
            tempVect = tempVect(isfinite(tempVect));

            s_fieldname = ['mean_',fields{i}];
            fd_summ.(s_fieldname) = mean(tempVect,'omitnan');

            s_fieldname = ['std_',fields{i}];
            fd_summ.(s_fieldname) = std(tempVect,'omitnan');
        end
    end
    fd_summ.metrics = metrics;
end