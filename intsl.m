function a=intsl(f) %实现left spectral integration matrix
k=size(f);
k=max(k);
a=zeros(1,k+1);
a(2)=(2*f(1)-f(3))/2;
for i=3:k-1
    a(i)=(f(i-1)-f(i+1))/(2*(i-1));
end
a(k)=f(k-1)/(2*(k-1));
a(k+1)=f(k)/(2*k);
for i=2:k+1
    a(1)=a(1)+(-1)^i*a(i);
end
a=a(1:k);