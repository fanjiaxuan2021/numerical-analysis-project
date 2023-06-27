function [matrix] = P(a,b,k,m)
t=1:m+1;
h=(b-a)/m;
z=a+h.*(t-1);

x=zeros(m,k);
for l=1:m
    x(l,:)=cheb(z(l),z(l+1),k);
end

w=zeros(m,k);
for u=1:m
    for v=1:k
        w(u,:)=weight(z(u),z(u+1),k);
    end
end

matrix1=zeros(m*k,m*k);
for i=1:m
    for j=1:k
        for r=1:i-1
            for s=1:k
                matrix1(k*(i-1)+j,k*(r-1)+s)=phil(x(i,j))*gl(x(r,s))*w(r,s);
            end
        end
        Q=Isl(k,z(i),z(i+1));
        for s=1:k
            matrix1(k*(i-1)+j,k*(i-1)+s)=phil(x(i,j))*gl(x(i,s))*Q(j,s);
        end
    end
end

matrix2=zeros(m*k,m*k);
for i=1:m
    for j=1:k
        R=Isr(k,z(i),z(i+1));
        for s=1:k
            matrix2(k*(i-1)+j,k*(i-1)+s)=phir(x(i,j))*gr(x(i,s))*R(j,s);
        end
        for r=i+1:m
            for s=1:k
                matrix2(k*(i-1)+j,k*(r-1)+s)=phir(x(i,j))*gr(x(r,s))*w(r,s);
            end
        end
    end
end

matrix0=eye(m*k);

matrix=matrix0+matrix1+matrix2;

end