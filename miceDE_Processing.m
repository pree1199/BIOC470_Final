function mgoi = miceDE_Processing(mgoi,nametoid)
%MICEDE_PROCESSING Summary of this function goes here
%   Detailed explanation goes here
%go through already paired ensembl id and gene names and assign values in
%mgoi accordingly
for i=1:size(mgoi,1)
    for j=1:size(nametoid,1)
        if isequal(mgoi(i,1),nametoid{j,2})
            mgoi(i,2)=nametoid{j,1};
        end
    end
end
%for ids not covered in excel sheet, manually enter them
for i=1:size(mgoi,1)
    if isempty(mgoi{i,2})
        mgoi(i,2)={input(['Find Ensembl Id of ' mgoi{i,1}])};
    end
end
%process table; remove any genes with empty IDs
toremove=zeros(size(mgoi,1),1);
size(toremove)
counter=1;
for i=1:size(mgoi,1)
    if isempty(mgoi{i,2})
        toremove(counter)=i;
        counter=counter+1;
    end
end
toremove = toremove(toremove ~= 0);
for i=1:size(toremove,1)
    mgoi(toremove(i),:) = [];
end
%if first 7 characters of ID is not 'ENSMUSG', set entry to 0x0 char
for i=1:size(mgoi,1)
    if ~isequal({'ENSMUSG'},extractBetween(mgoi{i,2},1,7))
        mgoi{i,2}='';
    end
end
%manually enter missing IDs
for i=1:size(mgoi,1)
    if isempty(mgoi{i,2})
        mgoi(i,2)={input(['Find Ensembl Id of ' mgoi{i,1}])};
    end
end
%if any IDs are still missing, remove the genes from mgoi
toremove=zeros(size(mgoi,1),1);
size(toremove)
counter=1;
for i=1:size(mgoi,1)
    if isempty(mgoi{i,2})
        toremove(counter)=i;
        counter=counter+1;
    end
end
toremove = toremove(toremove ~= 0);
for i=1:size(toremove,1)
    mgoi(toremove(i),:) = [];
end
%write mgoi into a text file to save information
writecell(mgoi,'mousegenes.txt')
end