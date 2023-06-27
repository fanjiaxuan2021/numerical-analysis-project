function[y]=G0(x,t)
if (x<t)
    y=gl(x)*gr(t)/2;
else
    y=gl(t)*gr(x)/2;
end
end