function subclust_mat = createSubClusterMatrix(clust_mat,subclustnum,pcaID)
%CREATESUBCLUSTERMATRIX Summary of this function goes here
%   Detailed explanation goes here
ref=unique(pcaID);
ref = ref(ref~=0);
idoi=ref(subclustnum);
indexes=find(pcaID==idoi);
subclust_mat=zeros(size(clust_mat,1),length(indexes));
for i=1:length(indexes)
    subclust_mat(:,i)=clust_mat(:,indexes(i));
end
end

