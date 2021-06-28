function A = arr4mat(f,p) %f is a multicomplex number
arr = f.zn;
for l = 1:p-1
    for i = 1: length(f.zn)
        for j = 1 : length(f.zn)
            b = j;
            a = i;
            x = arr;
            n = length(f.zn);
            while length(x) > 1
                if a <= n/2 & b <= n/2 % upper left
                    x = x(1:n/2);
                elseif a <= n/2 & b > n/2 % upper right
                    x = -x(n/2 +1 :end);
                    b = b - n/2;
                elseif a > n/2 & b <= n/2 % lower left
                    x = x(n/2 +1 :end);
                    a = a - n/2;
                elseif a > n/2 & b > n/2 % lower right
                    x = x(1:n/2);
                    a = a - n/2;
                    b = b - n/2;
                end
                n = n/2;
            end
            s(j)=x; % up to this point
        end 
        y(i) = dot(s,f.zn);
    end
    arr = y;
end
A = y;
end
