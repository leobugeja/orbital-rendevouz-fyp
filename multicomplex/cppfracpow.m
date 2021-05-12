function A = cppfracpow(f,p)

n = length(f);

de = f(1)*n;
nu = n*f(1)^2;
k = nu/de;

C=((1/k)*f)-inputconverter(1,[1:log2(n)],0);

A = inputconverter(1,[1:log2(n)],0) + p * C;

temp = p*(p-1);
s = zeros(1,length(C));
y = zeros(1,length(C));
for i = 2:15 
    
    arr = C;
    for pow = 1:i-1
        for m = 1: length(C)
            for j = 1 : length(C)
                b = j;
                a = m;
                x = arr;
                l = length(C);
                while length(x) > 1
                    if a <= l/2 & b <= l/2 % upper left
                        x = x(1:l/2);
                    elseif a <= l/2 & b > l/2 % upper right
                        x = -x(l/2 +1 :end);
                        b = b - l/2;
                    elseif a > l/2 & b <= l/2 % lower left
                        x = x(l/2 +1 :end);
                        a = a - l/2;
                    elseif a > l/2 & b > l/2 % lower right
                        x = x(1:l/2);
                        a = a - l/2;
                        b = b - l/2;
                    end
                    l = l/2;
                end
                s(j)=x; % up to this point
            end          
            y(m) = dot(s,C);
        end
        arr = y;
    end
    C1 = y;
    
    A = A + (temp/factorial(i)) * C1;
    
    temp = temp * (p - i);
end

A = multicomplex((k^p)*A);
end