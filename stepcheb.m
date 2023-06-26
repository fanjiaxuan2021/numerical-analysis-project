function F = stepcheb(b,u,k,x)
F=0;
m=length(b)-1;
for i=1:m
    if x>=b(i)&&x<b(i+1)
        z=zeros(1,k);
        for j=1:k
            z(j)=u(i*k-k+j);
        end
        y=chebf(k)*z;
        for j=1:k
            F=F+y(k)*chebfunction(k-1,2*(x-b(i))/(b(i+1)-b(i))-1);
        end
    end
end

