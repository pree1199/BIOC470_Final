function subclusterClassification(celltype,clust_sub,liver_table)
%SUBCLUSTERCLASSIFICATION %Using the average read counts table for each subcluster, we will now
%normalize the expression across the clusters to find which of the cell
%type marker genes specified in variable celltype are preferentially
%expressed in each cluster to ascertai cellular identity of each subcluster
%specify zmat to be an empty matrix such that in the end, each column will
%correspond to a subcluster and each row correspond to a gene in celltype
zmat=zeros(length(celltype),size(clust_sub,2));
counter=1;
for i=1:size(celltype,1)
    for j=1:size(liver_table,1)
        %find the proper row in liver_table corresponding to the gene
        %specified in celltype; then compute the zscore of the corresponding row in clust_sub 
        %to normalize the data and store it in zmat
        if isequal(liver_table{j,1},celltype(i,1))
            zmat(counter,:)=zscore(table2array(clust_sub(j,:)));
            counter=counter+1;
        end
    end
end
%create a heatmap based on our values in zmat
figure;
heatmap(clust_sub.Properties.VariableNames,celltype',zmat,'CellLabelColor','none')
c = parula;
colormap(c);
set(gca,'FontSize',8)
set(gcf,'color','w');
end

