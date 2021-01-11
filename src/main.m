%% Task variables
global num_m num_xi num_l num_g num_k num_M
num_m  = 2;
num_xi = 5;
num_l  = 10;
num_g  = 10;
num_k  = 2;
num_M  = 20;

%% Task matrices
syms m xi l g k M

A = [0, 1, 0, 0; ...
    (m*g - xi)/(m * l), 0, -(m*g - xi)/(m*l), 0; ...
    0, 0, 0, 1; ...
    0, 0, 0, -k/M];
    
B = [0; 0; 0; 1/M];


global mat_A mat_B
mat_A = double(subs(A, [m xi l g k M], [num_m num_xi num_l num_g num_k num_M]));
mat_B = double(subs(B, M, num_M));

%% Kalman criteria
Controlability = [B, A*B, A*A*B, A*A*A*B]

determinant = det(Controlability)

%% Characteristic polinom A + BK
syms lambda
K = sym('K', [1, 4]);
charpoly(A + B*K)

%%
tspan = [0, 100];
u = @(t, x) empty_control(t, x);
[ts, xs] = ode45(@(t, x)main_system(t, x, u), tspan, [0.01; 0; 3; 0]);

figure;
subplot(1, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\varphi$$', 'interpreter', 'latex');
subplot(1, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot\varphi$$', 'interpreter', 'latex');
%% Task 1. Full Feedback
%
%  D * K^T = d1 - d2, where d1 = poly(mu).
%
global vec_mu
vec_mu = [-2 + 1i, -2 - 1i, -3 + 0i, -3 + 0i];

d1 = poly(vec_mu);
d1 = d1(2:end).';
D = [0, 0, 0, -1/M;
    0, 0, -1/M, 0;
    0, (m*g - xi)/(M*m*l), 0, (m*g - xi)/(M*m*l);
    (m*g - xi)/(M*m*l), 0, (m*g - xi)/(M*m*l), 0];
d2 = [-k/M; (xi - m*g)/(m*l); k*(xi - m*g/(M*m*l)); 0];

mat_D = double(subs(D, [m xi l g k M], [num_m num_xi num_l num_g num_k num_M]));
vec_d = double(subs(d1 - d2, [m xi l g k M], [num_m num_xi num_l num_g num_k num_M]));

mat_K = linsolve(mat_D, vec_d).';

%% Linear Stabilizator with full feedback
tspan = [0, 10];
u = @(t, x) linear_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)linear_system(t, x, u), tspan, [1.1; 1.1; 0.1; 0.1]);
figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$l\varphi+s$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$l\dot\varphi+\dot s$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');
%%
tspan = [0, 10];
u = @(t, x) main_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)main_system(t, x, u), tspan, [0.1; 0.1; 0.1; 0.1]);
figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');

%% Task 2. Dynamic Feedback
global mat_C
mat_C = [1, 0, 0, 0; 0, 0, 1, 0];
Controlability = [mat_C.', A.' * mat_C.', A.' * A.' * mat_C.']


%%
q = [-1; 0; 0; 0];
A.' * q
A.' * A.' * q
A.' * A.' * A.' * q
%%
T = [-1 , 0,  -(m*g - xi)/(m*l),                 0;...
      0, -1,                  0, -(m*g - xi)/(m*l); ...
      0,  0,   (m*g - xi)/(m*l),                 0; ...
      0,  0,                  0,  (m*g - xi)/(m*l)];
T^(-1) * A.' * T

%%
mat_T = double(subs(T, [m xi l g k M], [num_m num_xi num_l num_g num_k num_M]));
p1 = num_k / num_M;
p2 = (num_xi - num_m * num_g)/(num_m * num_l);
p3 = num_k*(num_xi - num_m * num_g)/(num_M * num_m * num_l);
p4 = 0;

vec_nu = [-2, -1, -3, -4];

d1 = poly(vec_nu);
d1 = d1(2:end).';
mat_Gamma = [p1 - d1(1), p2 - d1(2), p3 - d1(3), p4 - d1(4); 0, 0, 0, 0];
mat_Key = [1 p1 p2 p3; 0 1 p1 p2; 0 0 1 p1; 0 0 0 1];
global mat_L
mat_L = (mat_Gamma*(mat_T * mat_Key)^(-1)).';
poly(mat_A.'-mat_C.'*mat_L.')
%%
tspan = [0, 8];
u = @(t, x) linear_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)observer_system(t, x, u), tspan, [0.11; 0.11; 0.1; 0.1; 1.12; 1.12; 0.11; 0.11]);

figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');

%%
figure;
subplot(2, 2, 1);
plot(ts, abs(xs(:, 1)-xs(:,5)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat z_1 - z_1|$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, abs(xs(:, 2)-xs(:,6)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat z_2 - z_2|$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, abs(xs(:, 3)-xs(:,7)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat z_3 - z_3|$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, abs(xs(:, 4)-xs(:,8)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat z_4 - z_4|$$', 'interpreter', 'latex');
%%
tspan = [0, 8];
linear_u = @(t, x) linear_control(t, x, mat_K);
main_u = @(t, x) main_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)main_observer_system(t, x, main_u, linear_u), tspan, [0.01; 0.01; 0.1; 0.1; num_l*0.01+0.01+0.01; num_l*0.01+0.01+0.01; 0.09; 0.09]);

%%
figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');
%%
figure;
subplot(2, 2, 1);
plot(ts, abs(xs(:, 1)-(xs(:,5) - xs(:,7))/num_l));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat x_1 - x_1|$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, abs(xs(:, 2)-(xs(:,6) - xs(:,8))/num_l));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat x_2 - x_2|$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, abs(xs(:, 3) - xs(:, 7)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat x_3 - x_3|$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, abs(xs(:, 4)-xs(:,8)));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$|\hat x_4 - x_4|$$', 'interpreter', 'latex');

%% Task 3. LQR Stabilizator
mat_Q = [1 0 0 0; 0 1 0 0; 0 0 100 0; 0 0 0 100];
mat_R = eye(1);

[mat_K, ~, ~] = lqr(mat_A, mat_B, mat_Q, mat_R, 0);
mat_K = - mat_K;

%%
tspan = [0, 50];
u = @(t, x) linear_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)linear_system(t, x, u), tspan, [1.1; 1.1; 0.1; 0.1]);
figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$l\varphi+s$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$l\dot\varphi+\dot s$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');
%%
tspan = [0, 3];
u = @(t, x) main_control(t, x, mat_K);
[ts, xs] = ode45(@(t, x)main_system(t, x, u), tspan, [1; 1; 0.1; 0.1]);

figure;
subplot(2, 2, 1);
plot(ts, xs(:, 1));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 2);
plot(ts, xs(:, 2));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot\varphi$$', 'interpreter', 'latex');
subplot(2, 2, 3);
plot(ts, xs(:, 3));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$s$$', 'interpreter', 'latex');
subplot(2, 2, 4);
plot(ts, xs(:, 4));
grid on;
xlabel('$$t$$', 'interpreter', 'latex');
ylabel('$$\dot s$$', 'interpreter', 'latex');

function u = empty_control(t, x)
    u = 0;
end
function u = linear_control(t, x, mat_K)
    u = mat_K * x;
end

function u = main_control(t, x, mat_K)
    global num_l
    
    z = [x(1)*num_l + x(3); x(2)*num_l + x(4); x(3); x(4)];
    u = mat_K * z;
end

function dxdt = linear_system(t, x, u)
    global mat_A mat_B
    
    dxdt = mat_A * x + mat_B * u(t, x);
end

function dxdt = main_system(t, x, u)
    global num_m num_xi num_l num_g num_k num_M
    dxdt = [x(2); ...
        -num_xi/(num_m*num_l)*x(1) + num_g/num_l*sin(x(1)) - 1/(num_M*num_l)*(u(t, x) - num_k * x(4))*cos(x(1)); ...
        x(4); ...
        u(t,x)/num_M - num_k/num_M*x(4)];
end

function dxdt = observer_system(t, x, u)
    global mat_A mat_B mat_C mat_L
    estimate = [x(1); x(2); x(3); x(4)]; % Real position
    observe = [x(5); x(6); x(7); x(8)];  % My position
    dxdt = [mat_A * estimate + mat_B * u(t,observe); ...
        mat_A * observe + mat_B * u(t,observe) + mat_L*(mat_C*estimate - mat_C*observe)];
end

function dxdt = main_observer_system(t, x, main_u, linear_u)
    global mat_A mat_B mat_C mat_L num_m num_xi num_l num_g num_k num_M num_l
    estimate = [x(1); x(2); x(3); x(4)]; % Real position
    imaging_estimate = [x(1)*num_l + x(3); x(2)*num_l + x(4); x(3); x(4)];
    observe = [x(5); x(6); x(7); x(8)];  % My position
    dxdt = [estimate(2); ...
        -num_xi/(num_m*num_l)*estimate(1) + num_g/num_l*sin(estimate(1)) - 1/(num_M*num_l)*(linear_u(t, observe) - num_k * x(4))*cos(estimate(1)); ...
        estimate(4); ...
        linear_u(t,observe)/num_M - num_k/num_M*estimate(4); ...
        mat_A * observe + mat_B * linear_u(t,observe) + mat_L*(mat_C*imaging_estimate - mat_C*observe)];
    
end