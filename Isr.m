function w=Isr(k,a,b)
%写出right spectral integration matrix的具体形式
w=zeros(k,k);
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
for i=1:k
    z=zeros(k,1);
    z(i)=1;
    [~,f]=cheint(z,k,a,b,t);
    c=intsr(f);
    for j=1:k
        w(j,i)=0;
        for l=1:k
            w(j,i)=w(j,i)+(b-a)/2*c(l)*cos((l-1)*acos(t(j)));
        end
    end
end