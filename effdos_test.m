
data = csvread('temp_dep_fit.csv');
temp = data(:,1);
pos = data(:,2);
zepos = pos(1,1);
shift = 0.*temp;
for i = 1:numel(pos)
    shift(i,1) = pos(i,1) - pos(1,1);
end



