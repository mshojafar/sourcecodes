function C = deletezero(A)

%delete columns of A containing all zero entries
set(0,'RecursionLimit',1000)
b = sum(abs(A));
if (sum(b == 0) >= 1)
    ind = find(b == 0,1,'first');
    A(:,ind) = [];
    C = deletezero(A);
else
    C = A;
end