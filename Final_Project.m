%read in data and perform geneCounts for general information
liver_table=readtable('GSE115469_Proc_Data.csv');
geneCounts(liver_table);
%%
%find patient 1-5 indices
p1se=[2 0];
p2se=[0 0];
p3se=[0 0];
p4se=[0 0];
p5se=[0 0];
column=2;
while isequal({'P1TLH'},extractBetween(liver_table.Properties.VariableNames{column},1,5))
    column=column+1;
    if column>size(liver_table,2)
        break
    end
end
p1se(2)=column-1;
p2se(1)=column;
while isequal({'P2TLH'},extractBetween(liver_table.Properties.VariableNames{column},1,5))
    column=column+1;
    if column>size(liver_table,2)
        break
    end
end
p2se(2)=column-1;
p3se(1)=column;
while isequal({'P3TLH'},extractBetween(liver_table.Properties.VariableNames{column},1,5))
    column=column+1;
    if column>size(liver_table,2)
        break
    end
end
p3se(2)=column-1;
p4se(1)=column;
while isequal({'P4TLH'},extractBetween(liver_table.Properties.VariableNames{column},1,5))
    column=column+1;
    if column>size(liver_table,2)
        break
    end
end
p4se(2)=column-1;
p5se(1)=column;
while isequal({'P5TLH'},extractBetween(liver_table.Properties.VariableNames{column},1,5))
    column=column+1;
    if column>size(liver_table,2)
        break
    end
end
p5se(2)=column-1;
%%
%perform PCA, deterine explained variance by num components, choose desired number of components
%perform PCA with desird number of components and then tsne on resulting
%score matrix
liver_mat=table2array(liver_table(:,2:size(liver_table,2)));%create corresponding matrix from table
liver_mat_transpose=liver_mat';
[coeff,sc,eig,tsquared,variance]=pca(liver_mat_transpose); %notice spca would be better
var_summed=length(variance); 
%the following code will calculate the cumulative variance explained by the
%first i components and plot the graph
for i=1:length(variance)
    var_summed(i)=sum(variance(1:i));
end
figure;
plot(1:length(variance),var_summed)
title('Variance Explained by Component');
xlabel('Number of Components');
ylabel('Explained Variance');
% we choose enough principal components such that 97% of the variance of
% the data is explained
numpcs=find(round(var_summed,2)==97.0000);
i=round(length(numpcs)/2);
num_pc=numpcs(i);
%perform pca by specifying the number of components
[coeff,sc,eig,tsquared,variance]=pca(liver_mat_transpose,'NumComponents',num_pc);
Y1 = tsne(sc); %perform a tsne visualization
figure;
%%
gscatter(Y1(:,1),Y1(:,2))
title('t-Distributed Stochastic Neighbor Embedding (tSNE) Visualization');
xlabel('tSNE_1');
ylabel('tSNE_2');
%plot tSNE color coded by source patient liver
plotByPatient(Y1,p1se,p2se,p3se,p4se,p5se)
%%
%based on tSNE output choose n=8,9 and perform kmeans clustering on the
%score matrix from PCA
numofclust = input('How many clusters does the tSNE plot suggest this data has?');
pcaID = kmeans(sc, numofclust);
% figure;
% histogram(pcaID) %use to visualize distribution of cells in each cluster
%occasionally the clustering will create "false" clusters of very few cells
%such that they most likley have no biological signficance; the following
%code will process our clusters such that we only consider biologically relevant
%clusters 
[N,~]=histcounts(pcaID);
bool=N>50; %if a cluster has fewe than 50 cells, we will deem it to be "false"
if length(find(bool))>(numofclust-2) %as we expect 8 or 9 clusters from our tSNE, if the number of true clusters is less than 8, we will run kmeans again; if we have at least 8 clusters, we will proceed with processing
% the following commented out lines can be used to visualize how our
% kmeans clusters line up with the tSNE clusters
%     figure;
%     hold on;
%     scatter(Y1(:,1),Y1(:,2),5,pcaID,'filled')
%     hold off;
    [clustered_table, newpcaID]=kmeans_Processing(pcaID,liver_mat);%clustered_table is average read count matrix; newpcaID sets the id of those that fail to meet our theshold to 0
    identifyCellTypes(clustered_table,liver_table);%create heatmap based on selected genes 
end
%%
%Visualize tSNE color coded by cluster
figure;
hold on;
gscatter(Y1(:,1),Y1(:,2),newpcaID)
[lgd,objd]=legend({'Unassigned','Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8'});
objhl = findobj(objd, 'type', 'line'); % objects of legend of type line
set(objhl, 'Markersize', 12);
hold off;
title('tSNE Visualization Colored By Cluster');
xlabel('tSNE_1');
ylabel('tSNE_2');
%%
%based on clustering results, we will find DE genes; to get a sense of what
%we are dealing with we execute the following code and obtain a graph documenting the number of 
% DE genes identified that are overxpressed and underexpressed according to
% our algorithm
%numgenes is a 8x2 matrix with each row corresponding to a cluster and row
%1 corresponds to DE genes that are overexpressed and row 2 corresponds
%with underexpressed DE genes; each entry is the number of genes
numgenes=zeros(size(clustered_table,2),2);
for i=1:size(clustered_table,2)
    [DE_up,DE_down] = findDEGenesBetween(clustered_table,i,5);
    numgenes(i,1)=length(DE_up);
    numgenes(i,2)=length(DE_down);
end
figure;
plot(1:size(clustered_table,2),numgenes(:,1),1:size(clustered_table,2),numgenes(:,2));
title('Number of Differentially Expressed (DE) Genes Identified per Cluster');
xlabel('Cluster');
ylabel('Number of DE Genes');
%results suggest our function is unable to properly find DE downregulated
%genes; this is most likely due to the great deal of 0's in our expression
%profiles of each gene due to the nature of a scRNA-seq experiment; we will
%have to ignore the DE_down genes for our analysis to be meaningful
%From here on out DE Genes will refer exclusively to genes that are
%expressed 5 fold highers in a cluster, compared to the other clusters
%We now write an excel file with the names of DE Genes for each cluster
%Find DE genes for first cluster and store it in a struct, DE_info
[DE_up,~] = findDEGenesBetween(clustered_table,1,5);
var_name{1}=cat(2,'Cluster_',num2str(1));        
DE_info=struct(var_name{1},{liver_table{DE_up,1}});
%Find DE genes for remaining clusters and store it in DE_info
for i=2:size(clustered_table,2)
    [DE_up,~] = findDEGenesBetween(clustered_table,i,5);
    var_name{1}=cat(2,'Cluster_',num2str(i));        
    DE_info.(var_name{1})=liver_table{DE_up,1};
end
%Write DE_info into an excel file
temp_DEtable=struct2table(DE_info,'AsArray',true);
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for i=1:size(temp_DEtable,2)
    temp_table=temp_DEtable(:,i);
    tempcellarr=temp_table{:,1};
    tempcellarr=cat(1,{cat(2,'Cluster_',num2str(i))},tempcellarr{1});
    writecell(tempcellarr, 'Differentially_Expressed_Genes.xlsx','Range',cat(2,Alphabet(i),num2str(1)))
end
%%
%Based on literature and cell identity results, we suspect some clusters
%might have additional subclusters or subpopulations worth investigating
%we color code our tSNE plot, this time using cells we know our
%differentially expresed according to certain cell types that our clusters
%have the identity of; the specific distribution patterns on the tsne
%suggest subpopulations do exist in these 3 clusters
%Cluster 1
tsneColorCodeGene('CALCRL',liver_table,Y1);
tsneColorCodeGene('VWF',liver_table,Y1);
%Cluster 7
tsneColorCodeGene('LYZ',liver_table,Y1);
tsneColorCodeGene('CD5L',liver_table,Y1);
%Cluster 8
tsneColorCodeGene('CD3D',liver_table,Y1);
tsneColorCodeGene('CD7',liver_table,Y1);
%%
%subcluster finding; after color coding by genes it appears that
%cluster1,7 and 8 contain subpopulations; we attempt to identify them here
%We find the subClusters of cluster 1 and store the average expression
%values of each subcluster in clust1_sub
[clust1_sub, ~, Ynew1]=findSubCluster(1,liver_mat,newpcaID);
%%celltype is our list of cell type marker genes of interest; we use this
%%in subclusterClassification
celltype={'MGP';'SPARCL1';'TM4SF1';'CLEC14A';'CCL14';
    'CLEC1B';'FCN2';'S100A13';'RAMP3';'INMT';'DNASE1L3';'LIFR';'KRT7';'KRT19';'SOX9';
    'EPCAM';'ACTA2';'COL1A1';'RBP1'};
subclusterClassification(celltype,clust1_sub,liver_table);
%now that we have ascertained the validity of our clustering, we attempt to
%find DE genes
%We read the identified DE genes for each cluster as a later input for
%newDEValid
DEgenesgen=readtable('Differentially_Expressed_Genes.xlsx');
%We find DE genes for subcluster1 and trim it using newDEValid; we store
%the information in a struct, DE_info1
[DE_up,~] = findDEGenesBetween(clust1_sub,1,5);
DE_up = newDEValid(DE_up,DEgenesgen,liver_table);
var_name{1}=cat(2,'Subcluster_',num2str(1));        
DE_info1=struct(var_name{1},{liver_table{DE_up,1}});
%We find DE genes for the remaining subclusters and trim it using newDEValid; we store
%the information in DE_info1
for i=2:size(clust1_sub,2)
    [DE_up,~] = findDEGenesBetween(clust1_sub,i,5);
    DE_up = newDEValid(DE_up,DEgenesgen,liver_table);
    var_name{1}=cat(2,'Subcluster_',num2str(i));        
    DE_info1.(var_name{1})=liver_table{DE_up,1};
end
%We write DE_info1 into a an excel sheet
temp_DEtable=struct2table(DE_info1,'AsArray',true);
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for i=1:size(temp_DEtable,2)
    temp_table=temp_DEtable(:,i);
    tempcellarr=temp_table{:,1};
    tempcellarr=cat(1,{cat(2,'Subcluster_',num2str(i))},tempcellarr{1});
    writecell(tempcellarr, 'DEGenes_Clust1.xlsx','Range',cat(2,Alphabet(i),num2str(1)))
end
%%
%We execute the same exact code as above when attempting to find
%subclusters for cluster 1 except we now use cluster 7
[clust7_sub, ~, Ynew7]=findSubCluster(7,liver_mat,newpcaID);
celltype={'LYZ';'S100A9';'HLA-DPB1';'CD5L';'MARCO';'VSIG4'};
subclusterClassification(celltype,clust7_sub,liver_table);
[DE_up,~] = findDEGenesBetween(clust7_sub,1,5);
DE_up = newDEValid(DE_up,DEgenesgen,liver_table,7);
var_name{1}=cat(2,'Subcluster_',num2str(1));        
DE_info7=struct(var_name{1},{liver_table{DE_up,1}});
for i=2:size(clust7_sub,2)
    [DE_up,~] = findDEGenesBetween(clust7_sub,i,5);
    DE_up = newDEValid(DE_up,DEgenesgen,liver_table,7);
    var_name{1}=cat(2,'Subcluster_',num2str(i));        
    DE_info7.(var_name{1})=liver_table{DE_up,1};
end
temp_DEtable=struct2table(DE_info7,'AsArray',true);
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for i=1:size(temp_DEtable,2)
    temp_table=temp_DEtable(:,i);
    tempcellarr=temp_table{:,1};
    tempcellarr=cat(1,{cat(2,'Subcluster_',num2str(i))},tempcellarr{1});
    writecell(tempcellarr, 'DEGenes_Clust7.xlsx','Range',cat(2,Alphabet(i),num2str(1)))
end
%%
%We execute the same exact code as above when attempting to find
%subclusters for cluster 1 except we now use cluster 8
[clust8_sub, ~, Ynew8]=findSubCluster(8,liver_mat,newpcaID);
celltype={'CD2';'CD3D';'TRAC';'GZMK';'GNLY';'PTGDS';'GZMB';'TRDC';'STMN1';'HMGB2';'TYMS';'MS4A1';'LTB';'CD37';'CD79B';'CD7';'KLRB1';'NKG7'};
subclusterClassification(celltype,clust8_sub,liver_table);
[DE_up,~] = findDEGenesBetween(clust8_sub,1,5);
DE_up = newDEValid(DE_up,DEgenesgen,liver_table,8);
var_name{1}=cat(2,'Subcluster_',num2str(1));        
DE_info8=struct(var_name{1},{liver_table{DE_up,1}});
for i=2:size(clust8_sub,2)
    [DE_up,~] = findDEGenesBetween(clust8_sub,i,5);
    DE_up = newDEValid(DE_up,DEgenesgen,liver_table,8);
    var_name{1}=cat(2,'Subcluster_',num2str(i));        
    DE_info8.(var_name{1})=liver_table{DE_up,1};
end
temp_DEtable=struct2table(DE_info8,'AsArray',true);
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for i=1:size(temp_DEtable,2)
    temp_table=temp_DEtable(:,i);
    tempcellarr=temp_table{:,1};
    tempcellarr=cat(1,{cat(2,'Subcluster_',num2str(i))},tempcellarr{1});
    writecell(tempcellarr, 'DEGenes_Clust8.xlsx','Range',cat(2,Alphabet(i),num2str(1)))
end
%%
%As Cluster 8 did not divide into the clusters I expected it to, I
%created a matrix composed of average read counts for DE genes and computed
%Pearson Correlation Coefficient; if the coeffcient is relatively large, it
%indicates that while genes are differentially expressed, the two
%populations found in cluster 8 are somehow related and may not be true sub
%clusters
lgoi=[DE_info8.Subcluster_1;DE_info8.Subcluster_2];
for i=1:length(lgoi)
    for j=1:size(liver_table,1)
    	if isequal(lgoi(i,1),liver_table{j,1})
            DEgene_clust8(i,1)=clust8_sub{j,1};
            DEgene_clust8(i,2)=clust8_sub{j,2};
        end
    end
end
corrcoef(DEgene_clust8) %r-value of 0.0017 suggests significant differences exist
%%
%Mice Human Correlation
mice_table=readtable('miceLayerCounts.xlsx');
numgenes=zeros(length(2:30),1);
for i=2:31
    mgoi=miceDE(mice_table,i); %each entry of mgoi is the name of a gene that is significantly DE based on our value for the threshold
    numgenes(i-1)=length(mgoi);
end
figure;
plot(2:31,numgenes(:,1));
title('Number of Differentially Expressed (DE) Genes Depending on p-Value Threshold');
xlabel('p-Value Threshold-Log Scale');
ylabel('Number of DE Genes');
%choose p=10^-30 as threshold; gives us 110 genes to work with
%each entry of mgoi is the name of a gene that is significantly DE
mgoi=miceDE(mice_table,30);
%we now have 110 genes of interest and seek to find human orthologs; z-score each
%one; R-score
nametoid=readtable('ensembleToId.xlsx');%nametoid will have Ensembl Id in column 1 and the gene name in column 2
mgoi=[mgoi repmat({''},size(mgoi,1),1)]; %second column of mgoi will have ensembl ids
%go through already paired ensembl id and gene names and assign values in
%mgoi accordingly using miceDE_Processing
mgoi=miceDE_Processing(mgoi,nametoid);
%find human orthologs using ensembl
%call the resulting output from ensembl 'mart_export_unique.txt'
orth_genes=readtable('mart_export_unique.txt');
%resize mgoi and add 2 columns to it, column 3 will correspond to name of
%the human ortholog gene, column 4 will be the Ensembl Id of the human gene
mgoi=[mgoi repmat({''},size(mgoi,1),2)];
for i=1:size(mgoi,1)
    for j=1:size(orth_genes,1)
        if isequal(mgoi(i,2),orth_genes{j,1})
            mgoi(i,3)=orth_genes{j,3};
            mgoi(i,4)=orth_genes{j,4};
            break
        end
    end
end
%if no human ortholog exists, that is if column 3 and 4 are empty for a row, remove the entry
%let the modified list be newmgoi
counter=0;
for i=1:size(mgoi,1)
    if ~isempty(mgoi{i,3})
        counter=counter+1;
    end
end
newmgoi=cell(counter,4);
counter=1;
for i=1:size(mgoi,1)
    if ~isempty(mgoi{i,3})
        newmgoi(counter,:)=mgoi(i,:);
        counter=counter+1;
    end
end
%check and see if human ortholog gene was detected in MacParland et al. sc-RNA seq
%experiment; discard those that were not; %let the modified list be newermgoi
counter=0;
for i=1:size(newmgoi,1)
    for j=1:size(liver_table,1)
        if isequal(liver_table{j,1},newmgoi(i,3))
            counter=counter+1;
            break;
        end
    end
end
newermgoi=cell(counter,4);
counter=1;
for i=1:size(newmgoi,1)
    for j=1:size(liver_table,1)
        if isequal(liver_table{j,1},newmgoi(i,3))
            newermgoi(counter,:)=newmgoi(i,:);
            counter=counter+1;
            break;
        end
    end
end
%create matrices where each column is a cell in the specified cluster and
%the rows are still genes
[~, newpcaID]=kmeans_Processing(pcaId5,liver_mat);
clust2_mat=createClusterMatrix(liver_mat,2,newpcaID);%matrix corresponding to cluster 2
clust3_mat=createClusterMatrix(liver_mat,3,newpcaID);%matrix corresponding to cluster 3
clust5_mat=createClusterMatrix(liver_mat,5,newpcaID);%matrix corresponding to cluster 5
%next steps isolate clust mat for the specified 85 genes in newermgoi; take z-score of each
%gene; cross correlate them
clust2mat_ortho = createClustMatOrthoGenes(clust2_mat,liver_table,newermgoi);%matrix corresponding to cluster 2 where the rows are only those corresponding to 85 genes in newermgoi, in the same order as well
clust3mat_ortho = createClustMatOrthoGenes(clust3_mat,liver_table,newermgoi);%matrix corresponding to cluster 3 where the rows are only those corresponding to 85 genes in newermgoi, in the same order as well
clust5mat_ortho = createClustMatOrthoGenes(clust5_mat,liver_table,newermgoi);%matrix corresponding to cluster5 where the rows are only those corresponding to 85 genes in newermgoi, in the same order as well
%take the mean of expression to get an average read count matrix for the
%specified 85 genes for each cluster
clust2_mean=mean(clust2mat_ortho,2);
clust3_mean=mean(clust3mat_ortho,2);
clust5_mean=mean(clust5mat_ortho,2);
%concatenate the 3 average read count matrices into one matrix
combined_mean_mat=[clust2_mean clust3_mean clust5_mean];
%take the z_score of the clusters n human for each gene and store in
%zmat_hum
zmat_hum=zscore(combined_mean_mat,1,2);
%take the z_score of the layers in mice for each gene and store in
%zmat_mouse
zmat_mouse=zeros(size(newermgoi,1),9);
for i=1:size(newermgoi,1)
    for j=3:size(mice_table,1)
        if isequal(mice_table{j,1},newermgoi(i,1))
             zmat_mouse(i,:)=zscore(str2double(mice_table{j,2:10}));
        end
    end
end
%perform the pearson correlation coefficient calculation and store values
%in r_mat; store the corresponding p-values in p_mat
r_mat=zeros(size(zmat_mouse,2),size(zmat_hum,2));
p_mat=zeros(size(zmat_mouse,2),size(zmat_hum,2));
for i=1:size(zmat_hum,2)
    for j=1:size(zmat_mouse,2)
        [temp1_mat,temp2_mat]=corrcoef(zmat_hum(:,i),zmat_mouse(:,j));
        r_mat(j,i)=temp1_mat(1,2);
        p_mat(j,i)=temp2_mat(1,2);    
    end
end
figure;
heatmap({'Cluster 2';'Cluster 3';'Cluster 5'},{'Layer 1';'Layer 2';'Layer 3';'Layer 4';'Layer 5';'Layer 6';'Layer 7';'Layer 8';'Layer 9'},r_mat);
colormap(winter)
sig=p_mat<0.05;