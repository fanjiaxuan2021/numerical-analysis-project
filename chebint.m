function [result] = chebint(y, a, b, k)
    w=weight(a,b,k);
    result=w*y';
end