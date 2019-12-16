function DE_up = newDEValid(DE_up,DEgenesgen,liver_table,clustnum)
%NEWDEVALID When identifying DE genes in subclusters for a particular cluster, we only seek genes
%that are DE expressed by the cluster when compared to other clusters;
%otehrwise too many non-biologically relevant genes come into play; this
%code ensures that DE genes identified by findDEGenesBetween are actually
%differentialy expressed by the cluster. The genes that are not are
%discarded
%store only genes that are DE expressed by the cluster in newDE_up
newDE_up=zeros(size(DE_up,1),size(DE_up,2));
%We iterate through the list of identified DE_up genes for a particular
%subcluster and for each we go through the list of DE genes for the
%cluster the subcluster is derived from
for i=1:length(DE_up)
    for j=1:size(DEgenesgen,1)
        %if the gene in DE_up is found in the list of differentially
        %expressed genes for the cluster (found in DEgenesgen{:,clustnum}), we will store it
        if isequal(liver_table{DE_up(i),1},DEgenesgen{j,clustnum})
            newDE_up(i)=DE_up(i);
            break;
        end
    end
end
%set DE_up to equal the nonzero entries of newDE_up
DE_up = newDE_up(newDE_up ~= 0);
end

