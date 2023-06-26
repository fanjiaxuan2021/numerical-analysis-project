function P=local(a,b,phil,phir,gl,gr,k)
P=zeros(k,k);
pl=eva(phil,a,b,k);
pr=eva(phir,a,b,k);
hl=eva(gl,a,b,k);
hr=eva(gr,a,b,k);
A=il(k);
B=ir(k);
for i=1:k
    for j=1:k
        P(i,j)=(b-a)*pl(j)*A(i,j)*hl(i)/2+(b-a)*pr(j)*B(i,j)*hr(i)/2;
    end
end
P=P+eye(k,k);
end

