function [expgene_counts,gnoi] = geneCounts(liver_table)
%GENECOUNTS Summary of this function goes here
%   Detailed explanation goes here
expgene_counts=zeros(1,size(liver_table,2)-1);
for i=2:size(liver_table,2)
    expgene_counts(1,i-1)=size(find(liver_table{:,i}),1)-sum(ismissing(liver_table{:,i}));
end
figure;
histogram(expgene_counts)
%mean is 1.3129*10^3
%std is 676.53
title('Number of Expressed Genes per Cell');
xlabel('Number of Expressed Genes');
ylabel('Number of Cells');
set(gcf,'color','w');
% gnoi=zeros(1,size(liver_table,1)-1);
% counter=1;
% for i=2:size(liver_table,1)
%     if isempty(find(liver_table{i,2:size(liver_table,2)},1))
%         gnoi(counter)=liver_table{i,1};
%         counter=counter+1;
%     end
% end
end

