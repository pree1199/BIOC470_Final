function clust_mat = createClusterMatrix(liver_mat,clustnum,newpcaID)
%We wish to create a matrix where each column is a cell corresponding to
%the cluster specified in clustnum
ref=unique(newpcaID);
ref = ref(ref~=0);
idoi=ref(clustnum);
indexes=find(newpcaID==idoi);
clust_mat=zeros(size(liver_mat,1),length(indexes));
for i=1:length(indexes)
    clust_mat(:,i)=liver_mat(:,indexes(i));
end
end

