function [clust_sub, trueclust, Ynew] = findSubCluster(clustnum,liver_mat,newpcaID)
%FINDSUBCLUSTER We perform tSNE on the data to gain a visulation and a
%guess as to how many subclusters exist. We prompt user input for number of
%clusters and perform k-means clustering with the input as the number of
%clusters
%We create a matrix containing the data only for cells in the desired
%cluster
cluster1_mat=createClusterMatrix(liver_mat,clustnum,newpcaID);
%Perform tSNE on this cluster (We do not perform PCA in order to maximize
%our ability to distinguish presumably similar, but distinct, subclusters,
%despite the additional computational expense)
Ynew = tsne(cluster1_mat');
figure;
gscatter(Ynew(:,1),Ynew(:,2))
title('tSNE Visualization');
xlabel('tSNE_1');
ylabel('tSNE_2');
%ask for user input as to how many clusters k-means should create
numofclust = input('How many clusters does the tSNE plot suggest this data has?');
numclust1=zeros(1,100);
%Run k-means enough times for a sensible clustering result to come about;
%we define sensibiliy as when each created cluster has at east 20 cells,
%suggesting some semblance of significance
for i=1:100
pcaIDclust1 = kmeans(cluster1_mat', numofclust);
[N,~]=histcounts(pcaIDclust1);
numclust1(i)=sum(N>20);
if numclust1(i)>=(numofclust)
    idoi=pcaIDclust1;
    break;
end
end
figure;
hold on;
gscatter(Ynew(:,1),Ynew(:,2),idoi)
title('tSNE Visualization Color Coded by Subcluster');
xlabel('tSNE_1');
ylabel('tSNE_2');
legend();
hold off;
%based on subclustering, create an average read counts table, called clust_sub, for each
%identified subcluster 
trueclust=nonzeros(unique(idoi));
subclust_mat=createSubClusterMatrix(cluster1_mat,1,idoi);
clust_sub=array2table(mean(subclust_mat,2),'VariableNames',{sprintf('Subcluster_%d',1)});
for i=2:length(trueclust)
    subclust_mat=createSubClusterMatrix(cluster1_mat,i,idoi);
    clust_sub=[clust_sub array2table(mean(subclust_mat,2),'VariableNames',{sprintf('Subcluster_%d',i)})];
end

end


