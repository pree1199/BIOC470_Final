function [DE_up,DE_down] = findDEGenesBetween(clustered_table,clustnum,foldchange)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%We seek to find DE genes that are over expressed (the index of which is stored in DE_up) and
%that are underexpressed (the index of which is stored in DE_down)
DE_up=zeros(size(clustered_table,1),1);
DE_down=zeros(size(clustered_table,1),1);
%we will iterate through each of the 20,007 genes in clustered_table; each
%row of clustered_table corresponds to a gene; we will store the index of
%the proper row in DE_up or DE_down if the gene is qualified
for i=1:size(clustered_table,1)
    %voi will be the 1x8 matrix of average expression of each cluster for
    %the current gene in question
    voi=clustered_table{i,:};
    %cval will be the expression value for the current cluster of interest
    cval=voi(clustnum);
    %remove the entry of the cluster of interest from voi
    voi(:,clustnum) = [];
    %if each entry of voi, when multiplied by the specified foldchange, is still less than or
    %equal to the average expression of the gene in the cluster of
    %interest, store it in DE_up
    voi_mod=voi*foldchange;
    change=cval>=voi_mod;
    if sum(change)==length(voi)
        DE_up(i)=i;
    else
        % if a gene is qualified to be part of DE_up, consider if it is
        % part of DE_down: 
        %if each entry of voi, when divided by the specified foldchange, is still greater than or
        %equal to the average expression of the gene in the cluster of
        %interest, store it in DE_down; disregard the genes for which the
        %cluster of interest has an expression of 0
        voi_mod=voi/foldchange;
        change=cval<=voi_mod;
        if ~cval==0
            if sum(change)==length(voi)
                DE_down(i)=i;
            end
        end
    end
end
DE_up=DE_up(DE_up ~= 0);
DE_down=DE_down(DE_down ~= 0);
end

