function A=chebb(k)
A=zeros(k,k);
for i=1:k
    for j=1:k
        A(i,j)=cos((i-1)*(2*k-2*j+1)*pi/2/k);
    end
end
end

