function[y]=weight(a,b,k)
y=(2/k)*ones(1,k);
q=zeros(1,k);
for i=1:floor((k-1)/2)
    for j=1:k
        q(j)=cos(i*(2*k-2*j+1)*pi/k);
    end
    q=q*(4/(k*(2*i-1)*(2*i+1)));
    y=y-q;
end   
y=y*((b-a)/2);
end