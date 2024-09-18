% find column idx
function col = findcol(spec, D2_Titles)
    col = find(ismember(D2_Titles,spec)==1);
end