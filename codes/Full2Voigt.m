function C=Full2Voigt(LAMBDA)
%% Documentation
% This is the script to convert the 3^4 matrix to Voigt notation.
% Added by Xiaoran Chen on 2018/09/20. 
% Matlab version 2016a.

%% Preallocate space for 6 by 6 tensor
C=zeros(6);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                m=i*kronecker(i,j)+(1-kronecker(i,j))*(9-i-j);
                n=k*kronecker(k,l)+(1-kronecker(k,l))*(9-k-l);
                C(m,n)=LAMBDA(i,j,k,l);
            end
        end
    end
end

