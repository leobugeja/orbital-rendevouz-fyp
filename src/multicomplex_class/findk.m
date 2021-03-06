%find k

function A = findk(m,p)


de=0;
for i = 1:length(m)
    
    de = de + m(i,i);
end

nu = sum(m(:).^2);

k = nu/de

I=eye(length(m));

C=(1/k)*m-I;

A = I + p*C

temp = p*(p-1);

for i = 2:3
    
    A = A + (temp/factorial(i)) * C^i;
    
    show = (temp/factorial(i)) * C^i
    
    temp = temp * (p - i);
    
end
k^p
A = (k^p)*A;
end
