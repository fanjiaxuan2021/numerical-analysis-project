function X = localinnerproduct(a,b,gl,gr,phil,phir,f,k)
X=zeros(1,6);
P=local(a,b,phil,phir,gl,gr,k);
pl=eva(phil,a,b,k);
pr=eva(phir,a,b,k);
hl=eva(gl,a,b,k);
hr=eva(gr,a,b,k);
f0=eva(f,a,b,k);
ql=P\pl;
qr=P\pr;
f1=P\f0;
y1=il(k)*(hl.*ql);
y2=il(k)*(hr.*ql);
y3=il(k)*(hl.*qr);
y4=il(k)*(hr.*qr);
y5=il(k)*(hl.*f1);
y6=il(k)*(hr.*f1);
X(1)=y1(k);X(2)=y2(k);X(3)=y3(k);X(4)=y4(k);X(5)=y5(k);X(6)=y6(k);
end

