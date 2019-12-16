function clustmat_ortho = createClustMatOrthoGenes(clustmat,liver_table,newermgoi)
%CREATECLUSTMATORTHOGENES We want to trim down clustmat to only include
%rows corresponding to genes of interest which is documented in newermgoi
%use newermgoi and livertable to find proper indexes and then apply index
%on clustmat to create clustmat_ortho
clustmat_ortho=zeros(size(newermgoi,1),size(clustmat,2));
counter=1;
for i=1:size(newermgoi,1)
    for j=1:size(liver_table,1)
        if isequal(liver_table{j,1},newermgoi(i,3))
            clustmat_ortho(counter,:)=clustmat(j,:);
            counter=counter+1;
            break;
        end
    end    
end
end

