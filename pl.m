a=-1;b=1;k=10;m=100;
matrix=P(a,b,k,m);
g=1:m+1;
h=(b-a)/m;
r=a+h.*(g-1);
x=zeros(m,k);
for l=1:m
    x(l,:)=cheb(r(l),r(l+1),k);
end
z=zeros(1,m*k);
for i=1:m
    for j=1:k
        z(k*(i-1)+j)=x(i,j);
    end
end
u=main(a,b,k,m);
plot(z,u)