function tau=cheb(a,b,k)
tau=zeros(k,1);
for j=1:k
    tau(j)=(b-a)*(1/2)*cos((2*k-2*j+1)*pi/2/k)+a/2+b/2;
end

