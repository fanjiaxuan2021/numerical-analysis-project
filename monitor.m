function S = monitor(a,b,phil,phir,gl,gr,f,k)
f0=eva(f,a,b,k);
P=local(a,b,phil,phir,gl,gr,k);
sigma=P\f0;
sigma0=chebf(k)*sigma;
S=abs(sigma0(k-1))+abs(sigma0(k)-sigma0(k-2));
end

