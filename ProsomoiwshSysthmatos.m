clc;

% Γνωστές παράμετροι για τα εκκρεμή
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
Lambda_1 = 1 / J1;
Lambda_2 = 1 / J2;
a1 = m1 * g * r;
a2 = -0.5 * k * r;
a3 = -0.5 * b * r;
a4 = m2 * g * r;
a5 = 0.5 * k * r;
a6 = 0.5 * b * r;
g11=1000;
g12=1000;
g13=1000;
g21=1000;
g22=1000;
g23=1000;


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

% Ορισμός διανυσμάτων θ1, θ2, φ1, φ2
theta_1 = [a1; a2; a3]; % Διάνυσμα 3x1 για θ1
theta_2 = [a4; a5; a6]; % Διάνυσμα 3x1 για θ2

% Φυσικές συναρτήσεις φ1 και φ2
phi_1 = @(y1, y2, y3, y4) [
    f1(y1);
    f2(y1, y3);
    f3(y1, y2, y3, y4)
];

phi_2 = @(y1, y2, y3, y4) [
    f4(y3);
    f5(y1, y3);
    f6(y1, y2, y3, y4)
];
% Παράμετροι για τα μοντέλα αναφοράς
Am = [0 1; -4 -2]; % Πίνακας Am
Bm = [0; 1];       % Πίνακας Bm

% Επιθυμητές τροχιές
theta1d = @(t) (pi / 6) * sin(2 * pi * t);
theta1d_dot = @(t) (pi / 3) * pi * cos(2 * pi * t);
theta2d = @(t) (pi/4)* sin(pi * t);
theta2d_dot = @(t) (pi / 4) * pi * cos(pi * t);

r1 = @(t) (2 * pi / 3) * (1 - pi^2) * sin(2 * pi * t) + ((2 * pi^2) / 3) * cos(2 * pi * t);
r2 = @(t) ((pi^2) / 2) * cos(pi * t) + (pi - ((pi^3) / 4)) * sin(pi * t);


u1 = @(y,t) (-[y(13); y(14)]' * [y(1); y(2)]) - [y(15); y(16); y(17)]' * phi_1(y(1), y(2), y(3), y(4))+ (y(18)*r1(t));
u2 = @(y,t) (-[y(19); y(20)]' * [y(3); y(4)]) - [y(21); y(22); y(23)]' * phi_2(y(1), y(2), y(3), y(4))+ (y(24)*r2(t));

P=[3/2 1/8 ; 1/8 5/16];
B=[0; 1];
% Αρχικές συνθήκες

X0_1 = [theta1d(0); theta1d_dot(0)];
X0_2 = [theta2d(0); theta2d_dot(0)];

K0 = [0; 0]; % Αρχική τιμή του K
THETA_hat0= [0; 0; 0]; % Αρχική τιμή του THETA_hat

% Συνδυασμός όλων των αρχικών συνθηκών
y0_combined = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; % [y1, y2, y3, y4, τ1, τ2, Xm1, Xm2, K]
sigma = 0.01;
% Χρονικό διάστημα προσομοίωσης
tspan = [0, 20];

% Διαφορικές εξισώσεις
dydt_combined = @(t, y) [
    % Δυναμική του πρώτου εκκρεμούς
    y(2);
    Lambda_1 * (theta_1' * phi_1(y(1), y(2), y(3), y(4)) - T1(y(5), y(6), y(2)) + u1(y,t) );
    y(4);
    Lambda_2 * (theta_2' * phi_2(y(1), y(2), y(3), y(4)) - T2(y(7),  y(8), y(4))+ u2(y,t));
    tau_1_dot(y(2), y(5)); % dτ1/dt
    y(6);
    tau_1_dot(y(4), y(7)); % dτ2/dt
    y(8);
    % Μοντέλο αναφοράς για το εκκρεμές 1
    Am * y(9:10) + Bm * r1(t);
    % Μοντέλο αναφοράς για το εκκρεμές 2
    Am * y(11:12) + Bm * r2(t);
    % Εξίσωση για την εξέλιξη του K
    g11*(([y(1) - y(9); y(2) - y(10)]' *P * B * [y(1); y(2)])- sigma * y(13:14));
    % Εξίσωση για την εξέλιξη του ΤΗΕΤΑ_hat
    g12*((phi_1(y(1), y(2), y(3), y(4)) * [y(1) - y(9); y(2) - y(10)]' *P *B)- sigma * y(15:17));
    % Εξίσωση για την εξέλιξη του L
    g13*(((-([y(1) - y(9); y(2) - y(10)]' *P * B * r1(t) ))- sigma * y(18))- sigma * y(18));
    % Εξίσωση για την εξέλιξη του K2
    g21*((([y(3) - y(11); y(4) - y(12)]' *P * B * [y(3); y(4)]))- sigma * y(19:20));
    % Εξίσωση για την εξέλιξη του ΤΗΕΤΑ_hat2
    g22*((phi_2(y(1), y(2), y(3), y(4)) * [y(3) - y(11); y(4) - y(12)]' *P *B)- sigma * y(21:23));
    % Εξίσωση για την εξέλιξη του L2
    g23*((-([y(3) - y(11); y(4) - y(12)]' *P * B * r2(t) ))- sigma * y(24));
    
];

% Επίλυση του συστήματος
[t, y] = ode45(dydt_combined, tspan, y0_combined);  

%% Διαγραμμα για θ1
theta_values = theta1d(t);% Υπολογισμός της συνάρτησης επιθυμητης τροχιας
theta1_values = y(:, 1);  % Εξαγωγή των τιμών για θ1 από τη λύση y(1)
theta1d_values = y(:,9);  % Υπολογισμός της επιθυμητής τροχιάς θd1(t) απο τη λυση y(9)

% Δημιουργία διαγράμματος
figure;
plot(t, theta1_values, 'b-', 'LineWidth', 2);
hold on;
plot(t, theta1d_values, 'r--', 'LineWidth', 2);  % κόκκινη γραμμή για θd1% μπλε γραμμή για θ1
hold on;
plot(t, theta_values, 'g-', 'LineWidth', 1.5);
xlabel('Χρόνος (s)');
ylabel('Γωνία (rad)');
legend('\theta_1', '\theta_{m1}','\theta_{d1}');
title('Γωνία θ1 και Επιθυμητή τροχιά θd1 κατά τη διάρκεια του χρόνου');
grid on;
hold off;

%% Διαγραμμα για θ1dot

theta1d_dot_values = theta1d_dot(t);% Υπολογισμός της συνάρτησης επιθυμητης τροχιας
theta1_dot_values = y(:, 2);  % Εξαγωγή των τιμών για θ1dot από τη λύση y(2)
theta1ddot_values = y(:,10);  % Υπολογισμός της επιθυμητής τροχιάς θd1(t)dot απο τη λυση y(10)

% Δημιουργία διαγράμματος
figure;
plot(t, theta1_dot_values, 'b-', 'LineWidth', 2);
hold on;
plot(t, theta1ddot_values, 'r--', 'LineWidth', 2);  
hold on;
plot(t, theta1d_dot_values, 'g-', 'LineWidth', 1.5);
xlabel('Χρόνος (s)');
ylabel('Γωνία (rad)');
legend('\theta_1 dot', '\theta_{m1} dot','\theta_{d1} dot');
title('Γωνία θ1dot και Επιθυμητή τροχιά θd1dot κατά τη διάρκεια του χρόνου');
grid on;
hold off;

%% Διαγραμμα για θ2
theta_2_values = theta2d(t);% Υπολογισμός της συνάρτησης επιθυμητης τροχιας
theta2_values = y(:, 3);  % Εξαγωγή των τιμών για θ1 από τη λύση y(1)
theta2d_values = y(:,11);  % Υπολογισμός της επιθυμητής τροχιάς θd1(t) απο τη λυση y(9)

% Δημιουργία διαγράμματος
figure;
plot(t, theta2_values, 'b-', 'LineWidth', 2);
hold on;
plot(t, theta2d_values, 'r--', 'LineWidth', 2);  
hold on;
plot(t, theta_2_values, 'g-', 'LineWidth', 1.5);
xlabel('Χρόνος (s)');
ylabel('Γωνία (rad)');
legend('\theta_2', '\theta_{m2}','\theta_{d2}');
title('Γωνία θ2 και Επιθυμητή τροχιά θd2 κατά τη διάρκεια του χρόνου');
grid on;
hold off;

%% Διαγραμμα για θ1dot

theta2d_dot_values = theta2d_dot(t);% Υπολογισμός της συνάρτησης επιθυμητης τροχιας
theta2_dot_values = y(:, 4);  % Εξαγωγή των τιμών για θ2dot από τη λύση y(4)
theta2ddot_values = y(:,12);  % Υπολογισμός της επιθυμητής τροχιάς θd2(t)dot απο τη λυση y(12)

% Δημιουργία διαγράμματος
figure;
plot(t, theta2_dot_values, 'b-', 'LineWidth', 2);
hold on;
plot(t, theta2ddot_values, 'r--', 'LineWidth', 2);  
hold on;
plot(t, theta2d_dot_values, 'g-', 'LineWidth', 1.5);
xlabel('Χρόνος (s)');
ylabel('Γωνία (rad)');
legend('\theta_2 dot', '\theta_{m2} dot','\theta_{d2} dot');
title('Γωνία θ2dot και Επιθυμητή τροχιά θd2dot κατά τη διάρκεια του χρόνου');
grid on;
hold off;

%% διαγραμματα σφαλματων 
% Υπολογισμός του σφάλματος e11
e11 = y(:, 1) - y(:, 9); % y(:, 1) είναι οι τιμές του y1 από τη λύση

% Σχεδιασμός της γραφικής παράστασης για e1
figure;
plot(t, e11, 'm-', 'LineWidth', 2); % μωβ γραμμή για το e11
xlabel('Χρόνος (s)');
ylabel('Σφάλμα e_11 (rad)');
title('Σφάλμα e_1 = y_1 - \theta_{1d} κατά τη διάρκεια του χρόνου');
grid on;

% Υπολογισμός του σφάλματος e21
e21 = y(:, 3) - y(:, 11); % y(:, 3) είναι οι τιμές του y3 από τη λύση

% Σχεδιασμός της γραφικής παράστασης για e21
figure;
plot(t, e21, 'c-', 'LineWidth', 2); % γαλάζια γραμμή για το e21
xlabel('Χρόνος (s)');
ylabel('Σφάλμα e_2 (rad)');
title('Σφάλμα e_2 = y_3 - \theta_{2d} κατά τη διάρκεια του χρόνου');
grid on;


% Υπολογισμός του σφάλματος e12
e12 = y(:, 2) - y(:, 10); % y(:, 2) είναι οι τιμές του y2 από τη λύση

% Σχεδιασμός της γραφικής παράστασης για e12
figure;
plot(t, e12, 'm-', 'LineWidth', 2); % μωβ γραμμή για το e12
xlabel('Χρόνος (s)');
ylabel('Σφάλμα e_12 (rad)');
title('Σφάλμα e_12 = y_2 - \theta_{1d} dot κατά τη διάρκεια του χρόνου');
grid on;

% Υπολογισμός του σφάλματος e22
e22 = y(:, 4) - y(:, 12); % y(:, 4) είναι οι τιμές του y4 από τη λύση

% Σχεδιασμός της γραφικής παράστασης για e22
figure;
plot(t, e22, 'c-', 'LineWidth', 2); % γαλάζια γραμμή για το e22
xlabel('Χρόνος (s)');
ylabel('Σφάλμα e_22 (rad)');
title('Σφάλμα e_22 = y_4 - \theta_{2d} dot κατά τη διάρκεια του χρόνου');
grid on;


%% Διαγραμματα ελεγκτων
% Υπολογισμός της u1(y, t) συναρτήσει του χρόνου
u1_values = zeros(length(t), 1); % Αρχικοποίηση πίνακα για τιμές u1
for i = 1:length(t)
    % Υπολογισμός u1 για κάθε χρονική στιγμή
    u1_values(i) = u1(y(i, :), t(i));
end

% Δημιουργία του figure
figure;
plot(t, u1_values, 'LineWidth', 2);
xlabel('Χρόνος t (s)');
ylabel('u_1(y, t)');
title('Η συνάρτηση u_1(y, t) συναρτήσει του χρόνου');
grid on;


% Υπολογισμός της u2(y, t) συναρτήσει του χρόνου
u2_values = zeros(length(t), 1); % Αρχικοποίηση πίνακα για τιμές u2
for i = 1:length(t)
    % Υπολογισμός u2 για κάθε χρονική στιγμή
    u2_values(i) = u2(y(i, :), t(i));
end

% Δημιουργία του figure για u2
figure;
plot(t, u2_values, 'LineWidth', 2);
xlabel('Χρόνος t (s)');
ylabel('u_2(y, t)');
title('Η συνάρτηση u_2(y, t) συναρτήσει του χρόνου');
grid on;