function [filtered_data] = neglect_outliers(data)
    filtered_data = data;
    [numsim, numf] = size(data);
    meansim = mean(data);
    stdsim = std(data);
    lowBound = meansim - 3*stdsim;
    upBound = meansim + 3*stdsim;
    for i = 1:numsim
        for j = 1:numf
            if data(i,j) < lowBound(j)
                filtered_data(i,j) = NaN;
            end
            if data(i,j) > upBound(j)
                filtered_data(i,j) = NaN;
            end 
        end
    end
end

