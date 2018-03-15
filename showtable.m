function f = showtable(tbl)

f = char(tbl{1}');
for k=2:length(tbl)
    f = [f char(tbl{k}')];
end
