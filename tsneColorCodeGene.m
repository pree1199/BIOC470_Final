function tsneColorCodeGene(goi,liver_table,Y1)
%TSNECOLORCODEGENE ColorCode our tSNE visualization based on each cell's
%expression on a specific gene, goi
%genenames will be the name of all genes tested for in our sc-RNA
genenames=liver_table{:,1};
for i=1:length(genenames)
    %when we find the row that corresponds to goi in liver_table we will
    %execute the if statement
    if isequal(goi,genenames{i})
        %we take each cell's expression value for this gene
        goi_table=liver_table(i,:);
        %we normalize expression using z-score
        norm_goi_mat=zscore(goi_table{1,2:size(goi_table,2)},1,2);
        %we create 4 categories based on normaization: bottom (signified by 1) corresponding
        %to the lowest expression, low (signified by 2), medium (signified by 3) and highest (signified by 4), corresponding
        %to the highest expression; based on where each cell falls in its
        %expression of goi, we assign the aforementioned values to visid where each entry of vsid corresponds
        %to each cell; the value in visid will be used in our gscatter
        %function when we visualize the tsne output
        cutoff=range(norm_goi_mat)/4;
        cutoffs=[min(norm_goi_mat),min(norm_goi_mat)+cutoff,min(norm_goi_mat)+2*cutoff,min(norm_goi_mat)+3*cutoff];
        visid=zeros(size(liver_table,2)-1,1);
        for j=2:size(liver_table,2)
            if norm_goi_mat(j-1)>=cutoffs(4)
                visid(j-1)=4;
            elseif norm_goi_mat(j-1)<cutoffs(4)&&norm_goi_mat(j-1)>=cutoffs(3)
                visid(j-1)=3;
            elseif norm_goi_mat(j-1)<cutoffs(3)&&norm_goi_mat(j-1)>=cutoffs(2)
                visid(j-1)=2;
            else
                visid(j-1)=1;
            end
        end
        break;
    end
end
figure;
cmap=colormap(cool(4));
gscatter(Y1(:,1),Y1(:,2),visid,cmap)
title(goi);
xlabel('tSNE_1');
ylabel('tSNE_2');
legend({'Bottom','Low','Medium','High'});
set(gcf,'color','w');
end

