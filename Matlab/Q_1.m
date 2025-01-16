clc;
clear;
close all;

%% First order limit (mean)
alpha = -1;
beta = 1;
h = 0.05;
N = 41;
x1 = alpha:0.01:beta;
x2 = alpha:0.01:beta;
[x1, x2] = meshgrid(x1, x2);
g_bar = zeros(N*N, 1);
e_i1 = zeros(N, 1);
e_i2 = zeros(N, 1);
num = 0;
den = 0;
k = 1;
trimf = @(x, abc) max(min((x - abc(1)) / (abc(2) - abc(1)), (abc(3) - x) / (abc(3) - abc(2))), 0);
for i1 = 2:N
    for i2 = 2:N
        e_i1(i1-1,1) = -1 + h*(i1-2);
        e_i2(i2-1,1) = -1 + h*(i2-2);
        if i1 == 2
            mu_A_x1 = trimf(x1, [-1, -1, -1+h]);
        elseif i1 == N
            mu_A_x1 = trimf(x1, [1-h, 1, 1]);
        else
            mu_A_x1 = trimf(x1, [-1+h*(i1-3), -1+h*(i1-2), -1+h*(i1-1)]);
        end
        if i2 == 2
            mu_A_x2 = trimf(x2, [-1, -1, -1+h]);
        elseif i2 == N
            mu_A_x2 = trimf(x2, [1-h, 1, 1]);
        else
            mu_A_x2 = trimf(x2, [-1+h*(i2-3), -1+h*(i2-2), -1+h*(i2-1)]);
        end
        g_bar(k,1) = 1 / (3 + e_i1(i1-1,1) + e_i2(i2-1,1));
        num = num + g_bar(k,1) * mu_A_x1 .* mu_A_x2;
        den = den + mu_A_x1 .* mu_A_x2;
        k = k + 1;
    end
end
f_x = num ./ den;
g_x = 1 ./ (3 + x1 + x2);
figure
surf(x1, x2, g_x, 'FaceColor', 'b', 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('g(x)', 'Interpreter', 'latex');
legend('$g(x)$', '$f(x)$', 'Interpreter', 'latex');
grid on;
figure
surf(x1, x2, f_x, 'FaceColor', 'r', 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('f(x)', 'Interpreter', 'latex');
legend('$f(x)$', 'Interpreter', 'latex');
grid on;
figure
surf(x1, x2, g_x - f_x, 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('Error', 'Interpreter', 'latex');
title('Error', 'Interpreter', 'latex');
colorbar;
grid on;

%% Second order limit (mean)
alpha = -1;
beta = 1;
h = 0.25;
N = 9;
x1 = alpha:0.01:beta;
x2 = alpha:0.01:beta;
[x1, x2] = meshgrid(x1, x2);
g_bar = zeros(N*N, 1);
e_i1 = zeros(N, 1);
e_i2 = zeros(N, 1);
num = 0;
den = 0;
k = 1;
% Define trimf function
trimf = @(x, abc) max(min((x - abc(1)) / (abc(2) - abc(1)), (abc(3) - x) / (abc(3) - abc(2))), 0);
% Loop to calculate memberships and g_bar
for i1 = 2:N
    for i2 = 2:N
        e_i1(i1-1,1) = -1 + h*(i1-2);
        e_i2(i2-1,1) = -1 + h*(i2-2);       
        if i1 == 2
            mu_A_x1 = trimf(x1, [-1, -1, -1+h]);
        elseif i1 == N
            mu_A_x1 = trimf(x1, [1-h, 1, 1]);
        else
            mu_A_x1 = trimf(x1, [-1+h*(i1-3), -1+h*(i1-2), -1+h*(i1-1)]);
        end       
        if i2 == 2
            mu_A_x2 = trimf(x2, [-1, -1, -1+h]);
        elseif i2 == N
            mu_A_x2 = trimf(x2, [1-h, 1, 1]);
        else
            mu_A_x2 = trimf(x2, [-1+h*(i2-3), -1+h*(i2-2), -1+h*(i2-1)]);
        end
        g_bar(k,1) = 1 / (3 + e_i1(i1-1,1) + e_i2(i2-1,1));
        num = num + g_bar(k,1) * mu_A_x1 .* mu_A_x2;
        den = den + mu_A_x1 .* mu_A_x2;
        k = k + 1;
    end
end
% Calculate f_x and g_x
f_x = num ./ den;
g_x = 1 ./ (3 + x1 + x2);
figure
surf(x1, x2, g_x, 'FaceColor', 'b', 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('g(x)', 'Interpreter', 'latex');
legend('$g(x)$', '$f(x)$', 'Interpreter', 'latex');
grid on;
figure
surf(x1, x2, f_x, 'FaceColor', 'r', 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('f(x)', 'Interpreter', 'latex');
legend('$f(x)$', 'Interpreter', 'latex');
grid on;
figure
surf(x1, x2, g_x - f_x, 'EdgeColor', 'none');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
zlabel('Error', 'Interpreter', 'latex');
title('Error', 'Interpreter', 'latex');
colorbar;
grid on;

%% First order limit (max)
alfa=-1;
beta=1;
h=0.05;
N=41;
x1=alfa:0.01:beta;
x2=x1;
[~,n1]=size(x1);
[~,n2]=size(x2);
e1=beta*ones(1,N+1);
e2=beta*ones(1,N+1);
for j=1:N
    e1(j)=alfa+h*(j-1);
    e2(j)=alfa+h*(j-1);
end
f_x=zeros(n1,n2);
for k1=1:n1
    for k2=1:n2

       i1=min(find(e1<=x1(1,k1),1,'last'),find(e1>=x1(1,k1),1));
       i2=min(find(e2<=x2(1,k2),1,'last'),find(e2>=x2(1,k2),1));
       
       if x1(1,k1)>=e1(1,i1)&& x1(1,k1)<=.5*(e1(1,i1)+e1(1,1+i1))&& x2(1,k2)>=e2(1,i2)&& x2(1,k2)<=.5*(e2(1,i2)+e2(1,1+i2))
           p=0;
           q=0;
       elseif x1(1,k1)>=e1(1,i1)&& x1(1,k1)<=.5*(e1(1,i1)+e1(1,1+i1))&& x2(1,k2)>=0.5*(e2(1,i2)+e2(1,1+i2))&& x2(1,k2)<=e2(1,1+i2)
           p=0;
           q=1; 
       elseif x1(1,k1)>=.5*(e1(1,i1)+e1(1,1+i1))&& x1(1,k1)<=e1(1,1+i1)&& x2(1,k2)>=e2(1,i2)&& x2(1,k2)<=0.5*(e2(1,i2)+e2(1,1+i2))
           p=1;
           q=0; 
       elseif x1(1,k1)>=.5*(e1(1,i1)+e1(1,1+i1))&& x1(1,k1)<=e1(1,1+i1)&& x2(1,k2)>=0.5*(e2(1,i2)+e2(1,1+i2))&& x2(1,k2)<=e2(1,1+i2)
           p=1;
           q=1; 
       end

       f_x(k1,k2)=1/(3+e1(1,i1+p)+e2(1,i2+q));
    end
end
[x1,x2]=meshgrid(x1,x2);
figure1 = figure('Color',[1 1 1]);
mesh(x1,x2,transpose(f_x));
xlabel('x_1')
ylabel('x_2')
zlabel('f(x)')

%% ُSecond order limit (max)
alfa=-1;
beta=1;
h=0.25;
N=9;
x1=alfa:0.01:beta;
x2=x1;
[~,n1]=size(x1);
[~,n2]=size(x2);
e1=beta*ones(1,N+1);
e2=beta*ones(1,N+1);
for j=1:N
    e1(j)=alfa+h*(j-1);
    e2(j)=alfa+h*(j-1);
end
f_x=zeros(n1,n2);
for k1=1:n1
    for k2=1:n2

       i1=min(find(e1<=x1(1,k1),1,'last'),find(e1>=x1(1,k1),1));
       i2=min(find(e2<=x2(1,k2),1,'last'),find(e2>=x2(1,k2),1));
       
       if x1(1,k1)>=e1(1,i1)&& x1(1,k1)<=.5*(e1(1,i1)+e1(1,1+i1))&& x2(1,k2)>=e2(1,i2)&& x2(1,k2)<=.5*(e2(1,i2)+e2(1,1+i2))
           p=0;
           q=0;
       elseif x1(1,k1)>=e1(1,i1)&& x1(1,k1)<=.5*(e1(1,i1)+e1(1,1+i1))&& x2(1,k2)>=0.5*(e2(1,i2)+e2(1,1+i2))&& x2(1,k2)<=e2(1,1+i2)
           p=0;
           q=1; 
       elseif x1(1,k1)>=.5*(e1(1,i1)+e1(1,1+i1))&& x1(1,k1)<=e1(1,1+i1)&& x2(1,k2)>=e2(1,i2)&& x2(1,k2)<=0.5*(e2(1,i2)+e2(1,1+i2))
           p=1;
           q=0; 
       elseif x1(1,k1)>=.5*(e1(1,i1)+e1(1,1+i1))&& x1(1,k1)<=e1(1,1+i1)&& x2(1,k2)>=0.5*(e2(1,i2)+e2(1,1+i2))&& x2(1,k2)<=e2(1,1+i2)
           p=1;
           q=1; 
       end

       f_x(k1,k2)=1/(3+e1(1,i1+p)+e2(1,i2+q));
    end
end
[x1,x2]=meshgrid(x1,x2);
figure1 = figure('Color',[1 1 1]);
mesh(x1,x2,transpose(f_x));
xlabel('x_1')
ylabel('x_2')
zlabel('f(x)')