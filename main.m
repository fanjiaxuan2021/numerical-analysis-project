function [result] = main(a,b,k,m)
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
y=f(z);

sigma=matrix\y';

% Compute J_l and J_r
J_l=zeros(1,m);
J_r=zeros(1,m);
J_l(1)=0;
J_r(m)=0;
for i=1:m-1
    J_l(i+1)=J_l(i)+chebint(gl(z(((i-1)*k+1):(i*k))).*sigma(((i-1)*k+1):(i*k))',r(i),r(i+1),k);
end
for i=m:-1:2
    J_r(i-1)=J_r(i)+chebint(gr(z(((i-1)*k+1):(i*k))).*sigma(((i-1)*k+1):(i*k))',r(i),r(i+1),k);
end

% Compute u
v=u(x);
for i=1:m
    M=Isl(k,r(i),r(i+1));
    N=Isr(k,r(i),r(i+1));
    for j=1:k
        v((i-1)*k+j)=v((i-1)*k+j)+gr(z((i-1)*k+j))/2*(J_l(i)+M(j,:)*(gl(z(((i-1)*k+1):(i*k)))'.*sigma(((i-1)*k+1):(i*k))))...
            +gl(z((i-1)*k+j))/2*(J_r(i)+N(j,:)*(gr(z(((i-1)*k+1):(i*k)))'.*sigma(((i-1)*k+1):(i*k))));
    end
end

result=zeros(1,m*k);
for i=1:m
    for j=1:k
        result(k*(i-1)+j)=v(i,j);
    end
end

end

