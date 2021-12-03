% littoral_cell_indices - figure out starting indices of cells
% write results in littoral_cells.csv
old_cell = 'blank'
cell_names = [];
cstart = [];
for i=1:length(transects)
    cell = transects(i).littoral_cell;
    if ~strcmp(cell,old_cell)
       cell_names = [cell_names; {cell}];
       cstart = [cstart; i];
    end
    old_cell = cell;
end
fid = fopen('littoral_cells.csv','w')
for i = 1:length(cell_names)
    fprintf(fid,'%d,%s\n',cstart(i),string(cell_names(i)));
end
fclose(fid)