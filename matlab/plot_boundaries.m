function [] = plot_boundaries( tbx, tby )
sizeT = size(tbx); sizeT = sizeT(1);
for i=1:sizeT
    bx = table2cell(tbx(i,:));
    bx = bx(~cellfun('isempty',bx));
    bx = cell2mat(cellfun(@str2num, bx, 'Uni',0));
    by = table2cell(tby(i,:));
    by = by(~cellfun('isempty',by));
    by = cell2mat(cellfun(@str2num, by, 'Uni',0));
    plot3(bx, by, zeros(size(bx)), '-b'); hold on;
end
end

