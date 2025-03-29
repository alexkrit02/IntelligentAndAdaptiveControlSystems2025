clc;
clear;

% Ορισμός παραμέτρων για το μοντέλο αναφοράς του εκκρεμούς 1
Am = [0 1; -4 -4]; % Πίνακας Am
Bm = [0; 1];       % Πίνακας Bm

% Επιθυμητή τροχιά θ1d και η παράγωγός της
theta1d = @(t) (pi/6) * sin(2*pi*t);
theta1d_dot = @(t) (pi/3) * pi * cos(2*pi*t);

% Είσοδος αναφοράς r1 για το πρώτο εκκρεμές
r1 = @(t) (4*pi/6)*(1-pi^2)*sin(2*pi*t) + ((4*pi^2)/6)*cos(2*pi*t);

% Αρχικές συνθήκες για το X_m1 = [θ1d \dot{θ1d}]^T
X0_1 = [theta1d(0); theta1d_dot(0)];

% Χρονικό διάστημα προσομοίωσης
tspan = [0 10];

% Σύστημα εξισώσεων για το μοντέλο αναφοράς του εκκρεμούς 1
dXmdt_1 = @(t, Xm) Am * Xm + Bm * r1(t);

% Επίλυση του συστήματος με τη μέθοδο ODE45 για το εκκρεμές 1
[t1, Xm1] = ode45(dXmdt_1, tspan, X0_1);

% Εξαγωγή αποτελεσμάτων για το εκκρεμές 1
theta1d_vals = Xm1(:, 1); % Τιμές θ1d
theta1d_dot_vals = Xm1(:, 2); % Τιμές \dot{θ1d}

% Ορισμός παραμέτρων για το μοντέλο αναφοράς του εκκρεμούς 2
% Επιθυμητή τροχιά θ2d και η παράγωγός της
r2 = @(t) ((pi^2)/2)*cos(pi*t) + (pi - (pi^3)/4)*sin(pi*t);

% Αρχικές συνθήκες για το X_m2 = [θ2d \dot{θ2d}]^T
X0_2 = [0; 0]; % Μηδενικές αρχικές συνθήκες για το δεύτερο εκκρεμές

% Σύστημα εξισώσεων για το μοντέλο αναφοράς του εκκρεμούς 2
dXmdt_2 = @(t, Xm) Am * Xm + Bm * r2(t);

% Επίλυση του συστήματος με τη μέθοδο ODE45 για το εκκρεμές 2
[t2, Xm2] = ode45(dXmdt_2, tspan, X0_2);

% Εξαγωγή αποτελεσμάτων για το εκκρεμές 2
theta2d_vals = Xm2(:, 1); % Τιμές θ2d
theta2d_dot_vals = Xm2(:, 2); % Τιμές \dot{θ2d}

% Οπτικοποίηση των αποτελεσμάτων για το εκκρεμές 1
figure;
subplot(2, 1, 1);
plot(t1, theta1d_vals, 'r', 'LineWidth', 1.5);
title('Θέση \theta_{1d} (Επιθυμητή)');
xlabel('Χρόνος (s)');
ylabel('\theta_{1d} (rad)');
grid on;

subplot(2, 1, 2);
plot(t1, theta1d_dot_vals, 'r--', 'LineWidth', 1.5);
title('Γωνιακή Ταχύτητα \theta_{1d} dot');
xlabel('Χρόνος (s)');
ylabel('\theta_{1d} dot (rad/s)');
grid on;

% Οπτικοποίηση των αποτελεσμάτων για το εκκρεμές 2
figure;
subplot(2, 1, 1);
plot(t2, theta2d_vals, 'b', 'LineWidth', 1.5);
title('Θέση \theta_{2d} (Επιθυμητή)');
xlabel('Χρόνος (s)');
ylabel('\theta_{2d} (rad)');
grid on;

subplot(2, 1, 2);
plot(t2, theta2d_dot_vals, 'b--', 'LineWidth', 1.5);
title('Γωνιακή Ταχύτητα \theta_{2d} dot');
xlabel('Χρόνος (s)');
ylabel('\theta_{2d} dot (rad/s)');
grid on;
