clc;

% Γνωστές παράμετροι
J1 = 0.5;
J2 = 0.625;
m1 = 2; 
m2 = 2.5;
g = 9.81; 
r = 0.5; 
d = 0.5;
l = 0.5;
k = 150; 
b = 1;
sigma0 = 1;
sigma1 = 1;
sigma2 = 1;
theta_s_dot = 0.1;
Ts = 2;
Tc = 1;
Lambda_1=1/J1;
Lambda_2=1/J2;
a1 = m1 * g * r;
a2 = -0.5 * k * r;
a3 = -0.5 * b * r;
a4 = m2 * g * r;
a5 = 0.5 * k * r;
a6 = 0.5 * b * r;

x = @(y1, y3) sqrt(d^2 + d*r*(sin(y1)-sin(y3)) + (r^2)*(1 - cos(y3-y1))/2);
numerator = @(y1, y3) (r/2) * (cos(y3) - cos(y3)); % Αριθμητής
denominator = @(y1, y3) d + (r/2) * (sin(y3) - sin(y3)); % Παρονομαστής
theta = @(y1, y3) atan(numerator(y1, y3) / denominator(y1, y3)); % Υπολογισμός θ
numerator1 = @(y1, y2, y3, y4) (cos(y1)*y2 - cos(y3)*y4)*d*r + 0.5*(r^2)*sin(y3-y1)*(y4-y2);
denominator1 = @(y1, y3) 2*x(y1, y3);
xdot = @(y1, y2, y3, y4) numerator1(y1, y2, y3, y4) / denominator1(y1, y3);

% Define functions for friction
tau_1_dot = @(y2, y5) y2 - (sigma0 * abs(y2) / (Tc + (Ts - Tc) * exp(-abs(y2) / theta_s_dot))) * y5;
tau_2_dot = @(y4, y7) y4 - (sigma0 * abs(y4) / (Tc + (Ts - Tc) * exp(-abs(y4) / theta_s_dot))) * y7;


T1 = @(y5, y6, y2) sigma0 * y5 + sigma1 * y6 + sigma2 * y2;
T2 = @(y7, y8, y4) sigma0 * y7 + sigma1 * y8 + sigma2 * y4;

f1 = @(y1) sin(y1);                  
f2 = @(y1, y3) (x(y1, y3) - l) * cos(y1 - theta(y1, y3));                
f3 = @(y1, y2, y3, y4) xdot(y1, y2, y3, y4) * cos(y1 - theta(y1, y3));
f4 = @(y3) sin(y3);                  
f5 = @(y1, y3) (x(y1, y3) - l) * cos(y3 - theta(y1, y3));                
f6 = @(y1, y2, y3, y4) xdot(y1, y2, y3, y4) * cos(y3 - theta(y1, y3));

dydt = @(t, y) [
    y(2); % dy1/dt
    Lambda_1*( a1 * f1(y(1)) + a2 * f2(y(1), y(3)) + a3 * f3(y(1), y(2), y(3), y(4)) - T1(y(5), y(6), y(2)));
    y(4); % dy3/dt
    Lambda_2 * (a4 * f4(y(3)) + a5 * f5(y(1), y(3)) + a6 * f6(y(1), y(2), y(3), y(4)) - T2(y(7),  y(8), y(4)));
    tau_1_dot(y(2), y(5)); % dτ1/dt
    y(6);
    tau_1_dot(y(4), y(7)); % dτ2/dt
    y(8);
];

% Initial conditions
y0 = [0.001; 0; 0.001; 0; 0; 0; 0; 0]; % [y1, y2, y3, y4, τ1, τ2]

% Time span
tspan = [0, 10];

% Solve ODE
[t, y] = ode45(dydt, tspan, y0);

% Plot results
figure;
subplot(2, 1, 1);
plot(t, y(:, 1), 'r', t, y(:, 3), 'b');
title('State Variables y1 and y3');
xlabel('Time (s)');
ylabel('y1, y3');
legend('y1', 'y3');

subplot(2, 1, 2);
plot(t, y(:, 2), 'r--', t, y(:, 4), 'b--');
title('State Variables y2 and y4');
xlabel('Time (s)');
ylabel('y2, y4');
legend('y2', 'y4');


