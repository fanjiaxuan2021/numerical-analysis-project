function Z = upward(X,Y)
Z=zeros(1,6);
Z(1)=(1-Y(1))*(X(1)-X(3)*Y(2))/(1-X(3)*Y(2))+Y(1);
Z(2)=Y(2)*(1-X(4))*(1-X(1))/(1-X(3)*Y(2))+X(2);
Z(3)=X(3)*(1-Y(1))*(1-Y(4))/(1-X(3)*Y(2))+Y(3);
Z(4)=(1-X(4))*(Y(4)-X(3)*Y(2))/(1-X(3)*Y(2))+X(4);
Z(5)=(1-Y(1))*X(5)/(1-X(3)*Y(2))+Y(5)+(Y(1)-1)*X(1)*Y(6)/(1-X(3)*Y(2));
Z(6)=(1-X(4))*Y(6)/(1-X(3)*Y(2))+X(6)+(X(4)-1)*Y(4)*X(5)/(1-X(3)*Y(2));
end