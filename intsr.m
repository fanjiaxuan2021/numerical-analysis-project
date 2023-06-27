function b=intsr(f) %实现right spectral integration matrix
k=size(f);
k=max(k);
b=zeros(1,k+1);
b(2)=(f(3)-2*f(1))/2;
for i=3:k-1
    b(i)=(f(i+1)-f(i-1))/(2*(i-1));
end
b(k)=-f(k-1)/(2*(k-1));
b(k+1)=-f(k)/(2*k);
for i=2:k+1
    b(1)=b(1)-b(i);
end
b=b(1:k);