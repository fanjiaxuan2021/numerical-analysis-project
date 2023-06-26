function [S1,S2] = downward(X,Y,T)
S1=zeros(1,2);S2=zeros(1,2);
E=[1,Y(2);X(3),1];
F=[T(2)*(1-Y(4))-Y(6);T(1)*(1-X(1))-X(5)];
S=E\F;
S1(1)=T(1);S2(2)=T(2);
S1(2)=S(1);S2(1)=S(2);
end

