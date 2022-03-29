%nested for loops

for ii = 1:10 % number of rows
    for jj = 1:15 %number of columns
        A(ii,jj) = ii + jj;
    end
end
disp(A)