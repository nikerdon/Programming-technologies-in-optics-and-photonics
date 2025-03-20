clc;
clear all;
close all;

a = 0;
b = 10;

a1 = 3;
a2 = 2;

c0 = 2;
c1 = 0.5;

Y = dsolve('D2y+3*Dy+2*y=cos(t)','y(0)=2','Dy(0)=0.5')
ezplot(Y,[a b])
hold on

n = 100;

x = linspace(a, b, n);
y = linspace(a, b, n);
u = zeros(size(x));
phi = zeros(size(x));

h = (b-a)/n;

for i = 1:n
    k1 = 0;
    k2 = 0;
    k3 = 0;

    for j = 1:i-1
        K = a1 + a2*(x(i) - y(j));
        f = cos(y(j)) - c1*a1 - (c1*y(j) + c0)*a2;
        k1 = k1 + K*f;
        phi(i) = -k1*h + cos(x(i)) - c1*a1 - (c1*y(j) + c0)*a2;
        k2 = k2 + K*phi(j);
        phi(i) = -k2*h + cos(x(i)) - c1*a1 - (c1*y(j) + c0)*a2;
        k3 = k3 + phi(j)*(x(i)-y(j));
    end
    
    u(i) = c0 + c1*x(i) + k3*h;
end

scatter(x,u) 
legend('Exact','Numerical')