function e=eva(p,a,b,k)
z=cheb(a,b,k);
e=zeros(k,1);
e(k)=p(z(k));
end
