function A = fracpow(f,p)

n = length(f.zn);

de = f.zn(1)*n;
nu = n*f.zn(1)^2;
k = nu/de;

C=((1/k)*f)-1;

A = 1 + p * C;
A = A.zn;
temp = p*(p-1);

for i = 2:50
    
    A = A + (temp/factorial(i)) * arr4mat(C,i);
    
    temp = temp * (p - i);
end

A = multicomplex((k^p)*A);
end