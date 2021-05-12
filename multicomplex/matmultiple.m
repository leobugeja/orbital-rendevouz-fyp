function x = matmultiple(f,p) %

x = f.zn;
for l = 1: p-1
    
    a = 0; b = 2; large = 0;
    for i = 1: length(f.zn)
        if i == 1
            s = x;
            init = x;
            prev = x;
            
        elseif rem(i,2)==0 %if row number is even
            temp = prev;
            for j = 1: length(x)
                if rem(j,2)~=0 %if the number of the element is odd
                    s(j)=temp(j+1);
                else
                    s(j)=temp(j-1);
                end
            end
            
        else
            m = length(x)/2;
            
            while m >= 2 % how many element swap each time
                if rem(i-1,m)== 0 | rem(i-2,m)== 0 %even number of row
                    
                    n = m;
                    
                    if n > a
                        arr = init;
                    elseif n > b & n < a | n == a/2
                        arr = large;
                    elseif n < a & n < b | n < a
                        arr = templarge;
                        %if
                    end
                    
                    s = [arr(n+1:2*n),arr(1:n)];
                    
                    while length(s) < length(f.zn)
                        temp = [arr(length(s)+n+1:length(s)+2*n),arr(length(s)+1:length(s)+n)];
                        s = cat(2,s,temp);
                    end
                    
                    if n == a/2 | n > a
                        a = n;
                        b = 2;
                        large = s;
                        templarge = s;
                    elseif n > b | n > c
                        templarge = s;
                    end
                    
                    if n < a & n > b
                        b = n;
                    end
                    c = n;
                    break
                end
                m=m/2;
            end
        end
        prev = s;
        
        for j = i+1 : length(x)
            s(j) = -s(j);
        end
        
        tem = 0;
        
        for j = 1: length(x)
            tem = tem + s(j) * f.zn(j);
        end
        x(i) = tem;
    end
end



