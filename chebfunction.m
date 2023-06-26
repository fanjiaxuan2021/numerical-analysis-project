function T = chebfunction(k,x)
if k==0
    T=1;
elseif k==1
    T=x;
elseif k>=2
    T1=chebfunction(k-1,x);
    T2=chebfunction(k-2,x);
    T=2*x*T1-T2;
end
end

