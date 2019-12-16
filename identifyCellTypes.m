function clustered_table = identifyCellTypes(clustered_table,liver_table)
%Using the average read counts table for each cluster, we will now
%normalize the expression across the clusters to find which of the cell
%type marker genes specified in variable celltype are preferentially
%expressed in each cluster to ascertai cellular identity of each cluster

%celltype is our list of cell type marker genes of interest
celltype={'CYP3A7';'CYP2A6';'CYP2A7';'SCD';'HMGCS1';'ACSS2';'TM7SF2';'SEC16B';
    'SLBP';'RND3';'PCK1';'BCHE';'G6PC';'ALDH6A1';'RPP25L';'HSD11B1';'HAMP';
    'GHR';'HPR';'GSTA2';'AKR1C1';'MASP2';'MGP';'SPARCL1';'TM4SF1';'CLEC14A';'CCL14';
    'CLEC1B';'FCN2';'S100A13';'RAMP3';'INMT';'DNASE1L3';'LIFR';'KRT7';'KRT19';'SOX9';
    'EPCAM';'ACTA2';'COL1A1';'RBP1';'LYZ';'S100A9';'HLA-DPB1';'CD5L';'MARCO';'VSIG4';
    'CD2';'CD3D';'TRAC';'GZMK';'GNLY';'PTGDS';'GZMB';'TRDC';'STMN1';'HMGB2';'TYMS';'MS4A1';
    'LTB';'CD37';'CD79B';'IGLC2';'IGHG1';'IGKC';'CD7';'KLRB1';'NKG7';'HBB';'CA1';'ALAS2'};
%specify zmat to be an empty matrix such that in the end, each column will
%correspond to a cluster and each row correspond to a gene in celltype
zmat=zeros(length(celltype),size(clustered_table,2));
counter=1;
for i=1:size(celltype,1)
    for j=1:size(liver_table,1)
        %find the proper row in liver_table corresponding to the gene
        %specified in celltype; then compute the zscore of the corresponding row in clustered_table 
        %to normalize the data and store it in zmat
        if isequal(liver_table{j,1},celltype(i,1))
            zmat(counter,:)=zscore(table2array(clustered_table(j,:)));
            counter=counter+1;
        end
    end
end
%create a heatmap based on our values in zmat
figure;
heatmap(clustered_table.Properties.VariableNames,celltype',zmat)
c = parula;
colormap(c);
set(gca,'FontSize',8)
end

