function w=Isl(k,a,b)
%写出left spectral integration matrix的具体形式,其中该矩阵作用在函数在切比雪夫结点上的值
w=zeros(k,k);
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
for i=1:k
    z=zeros(k,1);
    z(i)=1;
    [~,f]=cheint(z,k,a,b,t);
    c=intsl(f);
    for j=1:k
        w(j,i)=0;
        for l=1:k
            w(j,i)=w(j,i)+(b-a)/2*c(l)*cos((l-1)*acos(t(j)));
        end
    end
end
