function h = multipleTtest(data)
    nDatasets = length(data);
    nTimepoints = size(data{1},2);
    pVals = nan(nTimepoints,1);
    
    
    for i=1:nTimepoints
        dataset1 = data{1};
        dataset2 = data{2};
        [~,p]=ttest2(dataset1(:,i), dataset2(:,i),'Vartype','unequal');
        pVals(i) = p;
    end
    h = fdr_bh(pVals);
end