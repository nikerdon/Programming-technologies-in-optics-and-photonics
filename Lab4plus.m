clc;
clear all;
eps = 5;
tol = 10^(-eps); % заданная точность

%A = [16 3 2; 3 5 1; 2 1 10]; % задание матрицы A
A = [8 3 2; 3 6 1; 2 1 7];
A = [11 3 2; 3 15 1; 2 1 12];
%A = [4 1 -1; 1 4 -1; -1 -1 4];
a = 0; % левая граница
b = 20; % правая граница

n = size(A,1);
sol = [];

[V,D,W] = eig(A);
w1 = vpa(V(:,1)/V(3,1),5);
w2 = vpa(V(:,2)/V(3,2),5);
w3 = vpa(V(:,3)/V(3,3),5);

w1/norm(w1)
w2/norm(w2)
w3/norm(w3)

syms f(x)
f(x) = det(A - x*eye(n));
df = diff(f,x);

%метод Ньютона для поиска корней характеристического уравнения
for i=a:0.5:b 
    x0 = i;
    while true
        if df(x0)==0
            break
        end
        x1 = round(x0 - f(x0)/df(x0),eps);
        
        if abs(x0 - x1) < tol
            break;
        end
        x0 = x1;
    end

    if ismember(round(x1,5),sol)==false
        sol = [sol round(x1,5)];
    end
end
sol

%поиск собственных векторов
for i=1:size(sol,2)
    a = A - sol(i)*eye(n)
    f = [0;0;0];
    
    D = diag(diag(a));
    A1 = tril(a,-1);
    A2 = triu(a,1);
    
    %метод Якоби
    xn1 = ones(n,1);
    
    for j=1:50
        xn1 = -inv(D)*A1*xn1 - inv(D)*A2*xn1 + inv(D)*f;
    end
    xn1 = xn1/norm(xn1);
    xn1 = round(xn1 / xn1(3,1),eps);
    if norm(a * xn1) <= 0.1
        xn1/norm(xn1)
    end

    %метод Зейделя
    xn2 = ones(n,1);
    n2 = 0;
    for j=1:50
        xn2 = inv(D+A1)*(f-A2*xn2);
    end
    xn2 = xn2/norm(xn2);
    xn2 = round(xn2 / xn2(3,1),eps);
    if norm(a * xn2) <= 0.1
        xn2/norm(xn2)
    end
    
    %метод SOR
    xn3 = ones(n,1);
    n3 = 0;
    w = 0.9;
    I = eye(size(a,1));
    
    for j=1:50
        xn3 = inv(I+w*inv(D)*A1)*(((1-w)*I-w*inv(D)*A2)*xn3+w*inv(D)*f);
    end
    xn3 = xn3/norm(xn3);
    xn3 = round(xn3 / xn3(3,1),eps);
    if norm(a * xn3) <= 0.1
        xn3/norm(xn3)
    end
end