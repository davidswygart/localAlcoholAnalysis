c = linearizeMat(control);
i = linearizeMat(inject);

c.treatment(:) = "control";
c.row = num2str(c.row)+"c";
i.treatment(:) = "inject";
i.row = num2str(i.row)+"i";

finalTable = [c;i];
finalTable = renamevars(finalTable,"col","bin");
finalTable = renamevars(finalTable,"row","clusterID");
finalTable = renamevars(finalTable,"val", "spk_z");

finalTable.bin = num2str(finalTable.bin)+"b";

writetable(finalTable,'alcoholInject.csv')

function outputTable = linearizeMat(mat)
linMat = mat(:);
[row,col] = ind2sub(size(mat), 1:length(linMat));

outputTable = table();
outputTable.row = row';
outputTable.col = col';
outputTable.val = linMat;
end
