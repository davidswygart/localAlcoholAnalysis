function saveCsvForR_repeatedMeasuresMixedModel(mats,names, fname)
    if length(names) < length(mats)
        error('Must enter 1 treatment name for each matrix')
    end

    t = cell(length(mats));
    for i = 1:length(mats)
        t{i} = linearizeMat(mats{i},names{i});
    end
    finalTable = vertcat(t{:});
    writetable(finalTable,fname)
end

function out = linearizeMat(mat,name)
    linMat = mat(:);
    [row,col] = ind2sub(size(mat), 1:length(linMat));
    
    out = table();
    out.clusterID = num2str(row') + string(name);
    out.bin = num2str(col') + "b"; %just adding text to force R to consider as string instead of continuous variable
    out.treatment(:) = {name};
    out.spk = linMat;
end
