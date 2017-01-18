s = tdfread('~/Desktop/sf.ndv', '\t');
s = struct2table(s);
s = s(:, 14:end);
s=table2array(s);



