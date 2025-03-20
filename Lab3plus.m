clc;
clear all;
close all;

a1 = 1;
a2 = 2;
w1 = 2;
w2 = 5;

%точный график функции
x0 = linspace(-1,1,100);
y0 = a1*cos(w1*x0)+a2*sin(w2*x0);
figure
plot(x0,y0)
hold on

%приближение Фурье
N = 100;
f = @(x1) a1*cos(w1*x1)+a2*sin(w2*x1);
sn =  @(x1) 0;
for i=0:N
    I1 = @(x1) f(x1).*cos(pi*i*x1);
    an = integral(I1,-1,1);

    I2 = @(x1) f(x1).*sin(pi*i*x1);
    bn = integral(I2,-1,1);
    
    if i==0
        sn = @(x1) an/2;
    else
        sn =  @(x1) sn(x1) + an*cos(pi*i*x1) + bn*sin(pi*i*x1);
    end
end
x1 = linspace(-1,1,100);
plot(x1,sn(x1))

%приближение Чебышёва
N = 32;
f = @(t) a1*cos(w1*cos(t))+a2*sin(w2*cos(t));
sn2 =  @(t) 0;
for i=0:N   
    if i==0
        I1 = @(t) f(t);
        a0 = 1/pi*integral(I1,0,pi);
        sn2 = @(t) a0;
    else
        I2 = @(t) f(t).*cos(i*t);
        an = 2/pi*integral(I2,0,pi);
        sn2 =  @(t) sn2(t) + an*cos(i*t);
    end
end
x1 = linspace(-1,1,100);
legend('f(x)','Fourier')
hold off
figure
plot(x1,abs(sn2(acos(x1))-y0))

%hold on
%plot(x0,y0)
legend('error','Chebyshev')