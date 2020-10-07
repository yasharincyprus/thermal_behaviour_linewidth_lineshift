function y = phi(x)
    line_feats = csvread('temp_dep_fit.csv');
    temp = line_feats(:,1);
    k = 0.695;
    T = 88/k;
    format long;
    y = x^3 / ((exp(x) - 1)*(sqrt((T/temp(1,1)^2)) + x));
end