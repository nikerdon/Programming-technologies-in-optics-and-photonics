clc;
clear all;

%a = [10 3 0; 3 15 1; 0 1 7];
a = [15 2 3 4; 2 11 3 4; 4 3 10 2; 4 3 2 15];
f = [4; 11; 5; 7];
%f = [2; 12; 5];
epsilon = 10^(-4);
%x = [-29/977; 748/977; 591/977];
x = inv(a)*f;

D = diag(diag(a));
A1 = tril(a,-1);
A2 = triu(a,1);

%метод Якоби
xn1 = [0; 0; 0; 0];
n1 = 0;
while norm(x-xn1)/norm(x) > epsilon
    n1 = n1+1;
    xn1 = -inv(D)*A1*xn1 - inv(D)*A2*xn1 + inv(D)*f;
end
n1

%метод Зейделя
xn2 = [0; 0; 0; 0];
n2 = 0;
while norm(x-xn2)/norm(x) > epsilon
    n2 = n2+1;
    xn2 = inv(D+A1)*(f-A2*xn2);
end
n2

%метод SOR
xn3 = [0; 0; 0; 0];
n3 = 0;
w = 0.9;
I = eye(size(a,1));

while norm(x-xn3)/norm(x) > epsilon
    n3 = n3+1;
    xn3 = inv(I+w*inv(D)*A1)*(((1-w)*I-w*inv(D)*A2)*xn3+w*inv(D)*f);
end
n3