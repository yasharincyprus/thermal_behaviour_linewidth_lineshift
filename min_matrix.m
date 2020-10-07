function [row,column,value] = min_matrix(A)
[y,in] = min(A);
[value,column] = min(y);
[~,row] = min(A(:,column));
end