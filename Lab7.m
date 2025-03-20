clc;
clear all;

N = 1000000;
output1 = 0;

a = -1;
b = 1;
c = -1;
d = 1;

A = 7;
B = 4;
C = 6;

f=@(x,y) A*x.^2+B*y.^2+C;

%Точное значение
ymin = @(x) -sqrt(1 - x.^2);
ymax = @(x) sqrt(1 - x.^2);
Exact = integral2(f,-1,1,ymin,ymax)

%Теорема о среднем
n = 0;
output1 = 0;
for i=1:N
   y=c+(d-c)*rand();
   x=a+(b-a)*rand();
   if x^2+y^2<=1
     output1=output1+f(x,y);
     n = n+1;
   end
end
Sred = pi*output1/n

%Ищем значение максимума функции в области
fun = @(x) -(A*x(1)^2 + B*x(2)^2 + C);

x0 = [0.1, 0.1];

A1 = [];
B1 = [];
Aeq = [];
beq = [];
lb = [-1, -1];
ub = [1, 1];
nonlcon = @(x) deal(x(1)^2 + x(2)^2 - 1, []);

[x, M] = fmincon(fun, x0, A1, B1, Aeq, beq, lb, ub, nonlcon);
M = -M;

%Метод Монте-Карло
n=0;
for i=1:N
   y=c+(d-c)*rand();
   x=a+(b-a)*rand();
   z=M*rand();
   if x^2+y^2<=1 && z<=f(x,y)
      n = n+1;
   end
end

Monte_Karlo = 4*M*n/N