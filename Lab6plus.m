clc;
clear all;
close all;

A = -14;
B = 0;
C = -10;
D = 0;
E = 0;
F = -5;

syms f(x,y) dfx(x,y) dfy(x,y) grad

f(x,y) = A*x^2/2 + B*x*y + C*y^2/2 - D*x - E*y + F;
dfx(x,y) = A*x+B*y-D;
dfy(x,y) = C*y+B*x-E;
grad(x,y) = [dfx(x,y); dfy(x,y)];

ddfx = A;
ddfxy = B;
ddfy = C;

x = linspace(-5,5);
y = linspace(-5,5);
[X,Y] = meshgrid(x,y);
Z = A.*power(X,2)/2 + B*X.*Y + C.*power(Y,2)/2 - D*X - E*Y + F;
contour(X,Y,Z)
hold on

x0 = 0;
y0 = 0;
scatter(x0,y0,'g','filled')

l = 0.1;
R = 5;
eps = 0.001;
n = 0;

while true
    x1 = x0 - l*dfx(x0,y0);
    y1 = y0 - l*dfy(x0,y0);

    if abs(sqrt(x1^2+y1^2)-sqrt(x0^2+y0^2)) < eps | x1^2+y1^2 >= R^2
        break
    end

    x0 = x1;
    y0 = y1;
    scatter(x0,y0,'r','filled')

    n = n+1;
end
double(x0)
double(y0)

x0 = 2;
y0 = -3;
n = 0;

H = [ddfx ddfxy; ddfxy ddfy];

while true
    r = [x0;y0] - inv(H)*grad(x0,y0);

    x1 = r(1,1);
    y1 = r(2,1);

    if abs(sqrt(x1^2+y1^2)-sqrt(x0^2+y0^2)) < eps | x1^2+y1^2 >= R^2
        break
    end

    x0 = x1;
    y0 = y1;
    scatter(x0,y0,'b','filled')

    n = n+1;
end
double(x0)
double(y0)

fun = @(x) A*x(1)^2/2 + B*x(1)*x(2) + C*x(2)^2/2 - D*x(1) - E*x(2) + F;
x0 = [2,-3];
x = fminsearch(fun,x0)
