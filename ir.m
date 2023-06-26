function A=ir(k)
I=zeros(k,k);
for j=3:k-1
    I(j,j-1)=-1/(2*(j-1));
    I(j,j+1)=1/(2*(j-1));
end
if k>2
    I(k,k-1)=-1/(2*(k-1));
    I(2,3)=1/2;
end
I(2,1)=-1;
I(1,1)=1;
I(1,2)=-1/4;
for i=3:k
    I(1,i)=-1/(i*(i-2));
end
A=chebb(k)*I*chebf(k);
end

