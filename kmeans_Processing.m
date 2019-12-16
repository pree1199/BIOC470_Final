function [clustered_table,pcaID] = kmeans_Processing(pcaID,liver_mat)
%We will create an average read counts table (clustered_table) and modify
%pcaID such that cells in a "false" cluster have their cluster id set to 0
r=range(pcaID); %subtracts max value of pcaID with min value; equivalent to find number of clusters
%the ith entry of instances will correspond to the number of cells in
%cluster i
instances=zeros(1,1+r);
for j=1:length(pcaID)
    for i=1:1+r
        if pcaID(j)==i
            instances(i)=instances(i)+1;
        end
    end
end
%the ith entry of instances_bool will equal 0 if cluster i has less than 50
%entries
instances_bool=instances>50;
%noi will be the indices of "false" clusters, or those that have less than
%50 entries
noi=find(~instances_bool);
%We now go through pcaID and set the cells assigned to a "false" cluster to
%0 so that they are now unassigned; later in our analysis we will assign them to a true
%cluster
for j=1:length(pcaID)
    for i=1:length(noi)
        if pcaID(j)==noi(i)
            pcaID(j)=0;
        end
    end
end
%initialize variable clustered_table
clustered_table=table;
%numclust will keep track of how many "true cluster" we have
numclust=0;
%i will be our cluster number
for i=1:1+r
    %creat a temporary matrix with each column corresponding to a cell in
    %cluster i
   temp_mat=zeros(size(liver_mat,1),instances(i));
   counter=1;
   %if the cluster is "true" we will  will find the cells in that cluster
   %and store their value for each gene in temp_mat
   if instances_bool(i)
       numclust=numclust+1;
        for j=1:length(pcaID)
            if pcaID(j)==i
                temp_mat(:,counter)=liver_mat(:,j);
                counter=counter+1;
            end
        end
        %take the mean of temp_mat to find the average expression value of
        %the cluster for each gene; convert the average matrix to a table
        %and concatenate it onto clustered_table with the proper Variable
        %Name
        temp_table=array2table(mean(temp_mat,2));
        temp_table.Properties.VariableNames{'Var1'} = sprintf('Cluster_%d', numclust);
        clustered_table=[clustered_table temp_table];
   end
end
end

