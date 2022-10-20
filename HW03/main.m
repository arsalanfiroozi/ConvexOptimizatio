%% Q3 a)
clearvars
clc
delta1 = 0;
delta2 = 0;
n = 2;
P = [1 -0.5;-0.5 2]; 
cvx_begin
    variable x(n)
    dual variable y1
    dual variable y2
    dual variable y3
    minimize( quad_form(x,P) - x(1) )
    subject to
        y1: x(1) + 2*x(2) <= -2 + delta1
        y2: x(1) - 4*x(2) <= -3 + delta2
        y3: 5*x(1) + 76*x(2) <= 1
cvx_end
quad_form(x,P) - x(1)
x
[y1 y2 y3]
x(1) + 2*x(2)
x(1) - 4*x(2)
5*x(1) + 76*x(2)
y1*(x(1) + 2*x(2)+2)
y2*(x(1) - 4*x(2)+3)
y3*(5*x(1) + 76*x(2)-1)
2*x(1) - x(2) - 1 + (y1+y2+5*y3)
4*x(2)-x(1)+(2*y1-4*y2+76*y3)
% Delta 0 0 8.22 -2.3 0.1667
%       0 -0.1 8.7 -2.4 0.17
%       0 0.1 7.98 -2.3 0.15
%       -0.1 0 8.57 -2.4 0.15
%       -0.1 -0.1 8.82 -2.4 0.1667
%       -0.1 +0.1 8.32 -2.3667 0.1333
%       +0.1 0 8.22 -2.3 0.1667
%       +0.1 -0.1 8.71 -2.4 0.1719
%       +0.1 +0.1 7.75 -2.25 0.1615
%% Q3 b)
clearvars
clc
deltas1 = [0 -0.1 0.1];
deltas2 = [0 -0.1 0.1];
n = 2;
P = [1 -0.5;...
     -0.5 2]; 
Y = zeros(9,3); 
X = zeros(9,2);
R = zeros(9,2);
k = 1;
cvx_begin
    variable x(n)
    dual variable y1
    dual variable y2
    dual variable y3
    minimize( quad_form(x,P) - x(1) )
    subject to
        y1: x(1) + 2*x(2) <= -2
        y2: x(1) - 4*x(2) <= -3
        y3: 5*x(1) + 76*x(2) <= 1
cvx_end
ps = quad_form(x,P) - x(1) ;
for delta1 = deltas1 
    for delta2 = deltas2
        [delta1 delta2]
        cvx_begin
            variable x(n)
            minimize( quad_form(x,P) - x(1) )
            subject to
                x(1) + 2*x(2) <= -2 + delta1
                x(1) - 4*x(2) <= -3 + delta2
                5*x(1) + 76*x(2) <= 1
        cvx_end
        Y(k,:) = [y1 y2 y3];
        X(k,:) = x;
        R(k,1) = quad_form(x,P) - x(1);
        R(k,2) = ps - y1 * delta1 - y2 * delta2;
        k = k + 1;
    end
end
%% Q4 b)
clearvars
clc
cvx_begin
    variable x(n)
    dual variable y1
    minimize( exp(-x(1)) )
    subject to
        y1: x(1)^2 <= 0
%         y1: x(1)^2/x(2) <= 0 % This line makes error
        x(2) > 0
cvx_end
x
y1
%% Q5_1
clearvars
clc
p = 6; n = 9; eps = .1;

% Incidence matrix
A = [-1 +1  0  0  0  0  0  0 -1;
     +1  0 -1  0  0  0  0 -1  0;
      0 +1 +1 -1  0  0  0  0  0;
      0  0  0 +1 -1 -1  0  0  0;
      0  0  0  0 +1  0 +1  0 +1;
      0  0  0  0  0 +1 -1 +1  0];

% Source rates
s = [1;0;0;0;-1;0];
t = [0;1;0;0;0;-1];

cvx_begin
    variables x(n) y(n)
    minimize( phi(x,y) )
    subject to
        A * x + s == 0
        A * y + t == 0
        x > 0
        y > 0
cvx_end
phi(x,y)
% 1.8466
% x: [[0.131183848178559,5.75902202881727e-13,0.130020273513626,0.130020273514049,0.130020273362228,1.51558726248142e-10,0.00116357481610756,0.00116357466486557,0.868816151821982]]
% y: [[5.84365514477061e-13,1.20911515194703e-12,0.191164326950314,0.191164326951210,2.09247032906923e-13,0.191164326950620,2.08531835940276e-13,0.808835673050015,1.23594074947928e-13]]
%% Q5_2
clearvars
addpath(genpath('matlabGiftiCifti'))
clc
p = 6; n = 9; eps = .1;

% Incidence matrix
A = [-1 +1  0  0  0  0  0  0 -1;
     +1  0 -1  0  0  0  0 -1  0;
      0 +1 +1 -1  0  0  0  0  0;
      0  0  0 +1 -1 -1  0  0  0;
      0  0  0  0 +1  0 +1  0 +1;
      0  0  0  0  0 +1 -1 +1  0];

% Source rates
s = [1;0;0;0;-1;0];
t = [0;1;0;0;0;-1];

step_size = 0.1;
num_iterations = 200;

Y1 = zeros(1,p); Y2 = zeros(1,p);
x = zeros(n,1); y = zeros(n,1);
g0 = 0;
LB = [];
R1 = [];
R2 = [];
for i=1:num_iterations
    i
    for j=1:n
        [x1,y1] = quad2_min(eps,Y1*A(:,j),Y2*A(:,j));
        x(j) = x1;
        y(j) = y1;
        g(j) = phi(x1,y1) + Y1*A(:,j)*x1 + Y2*A(:,j)*y1;
    end
    LB = [LB phi(x,y) + Y1*A*x + Y2*A*y + Y1*s + Y2*t];
    R1 = [R1 norm(A*x+s)];
    R2 = [R2 norm(A*y+t)];
    g1 = phi(x,y) + Y1*A*x + Y2*A*y + Y1*s + Y2*t;
    Y1 = Y1 + transpose(step_size * (A*x + s));
    Y2 = Y2 + transpose(step_size * (A*y + t));
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(LB)
xlabel('Iterations');
ylabel('Lower Bound');
export_fig('1.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
semilogy(R1)
hold on
semilogy(R2)
box off
legend({ '||Ax+s||' , '||Ay+t||' });
xlabel('Iterations');
ylabel('Residual');
export_fig('2.png','-r600');
%% Q6 - a
clearvars
clc

step_size = 1;

Xs = [];
x = 1.1;
f = [];
Df = [];
for i = 1:10
    df = (exp(x) - exp(-x))/(exp(x) + exp(-x));
    Df = [Df df];
    f = [f log10(exp(x) + exp(-x))];
    d2f = 1+ 2* exp(-x)/(exp(x) + exp(-x)) + 2* exp(-x)/(exp(x) + exp(-x))^2 * (exp(x) - exp(-x));
    dx = - df / d2f;
    x = x + dx;
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(f)
hold on
plot(Df)
box off
legend({ 'f' , 'df' });
xlabel('Iterations');
export_fig('3_1.png','-r600');
%% Q6 - b
clearvars
clc

step_size = 1;

Xs = [];
x = 3;
f = [];
Df = [];
for i = 1:100
    df = 1-1/x;
    Df = [Df df];
    f = [f -log10(x)+x];
    d2f = 1/x^2;
    dx = - df / d2f;
    x = x + dx;
    if x<0
        x=0.0001;
    end
    Xs = [Xs x];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(f)
hold on
legend({ 'f' });
xlabel('Iterations');
export_fig('4_0.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(Df)
hold on
legend({ 'df' });
xlabel('Iterations');
export_fig('4_1.png','-r600');