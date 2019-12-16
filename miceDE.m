function mgoi = miceDE(mice_table,thresh)
%MICEDE Finds genes with a p-value less than 1*10^-thresh and returns it in
%a nx1 cell matrix where entry is the gene name
%let mgoi be a cell matrix with enough rows to accomodate all of the genes
%in mice_table
mgoi=cell(size(mice_table,1)-2,1);
counter=1;
%iterate 
for i=3:size(mice_table,1)
    tempvar=mice_table{i,20};
    if str2double(tempvar{1,1})<(1*10^(-thresh))
        mgoi(counter)=mice_table{i,1};
        counter=counter+1;
    end
end
mgoi=mgoi(~cellfun('isempty',mgoi));
end

