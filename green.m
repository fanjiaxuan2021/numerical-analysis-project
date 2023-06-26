function [gl,gr,g0l,g0r] = green(a,b,x00,x01,x10,x11)
if abs(x00)>=abs(x01)||abs(x10)>=abs(x11)
    gl=@(x) x00*(x-a)-x01;
    gr=@(x) x10*(x-b)-x11;
    g0l=@(x) x00;
    g0r=@(x) x10;
else
    gl=@(x) x01*cosh(x-a)-x00*sinh(x-a);
    gr=@(x) x11*cosh(x-b)-x10*sinh(x-b);
    g0l=@(x) x01*sinh(x-a)-x00*cosh(x-a);
    g0r=@(x) x11*sinh(x-b)-x10*cosh(x-b);
end

