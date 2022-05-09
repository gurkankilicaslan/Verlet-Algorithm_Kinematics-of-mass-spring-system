clc;
clear;
clear all;
dt = 0.1;
A = 1000;
x = zeros(3,A);
v = zeros(3,A);
a = zeros(3,A);

xmin=-2; xmax=2;
vmin=-2; vmax=2;
random_displacements = xmin + rand(1,3)*(xmax-xmin);
random_velocities = vmin + rand(1,3)*(vmax-vmin);

x(1,1) = random_displacements(1);
x(2,1) = random_displacements(2);
x(3,1) = random_displacements(3);

v(1,1) = random_velocities(1);
v(2,1) = random_velocities(2);
v(3,1) = random_velocities(3);

disp("These are the initial positions:")
disp(random_displacements)
disp("These are the initial velocities:")
disp(random_velocities)


m = 1;
k = 1;

k1 = k;
k2 = k;
m1 = m;
m2 = 3*m;
m3 = m;

%Analytical Solution using eigenvalues and eigenvectors
sm = sqrt([m1 m2 m3]');

Dynamic_Matrix = [
k1/sqrt(m1*m1) -k1/sqrt(m1*m2) 0
-k1/sqrt(m2*m1) (k1+k2)/sqrt(m2*m2) -k2/sqrt(m2*m3)
0 -k2/sqrt(m3*m2) k2/sqrt(m3*m3)
];

[V,D] = eig(Dynamic_Matrix);

w = sqrt(diag(D));
c_1 = V'*(sm.*x(:,1));
c_2 = (V'*(sm.*v(:,1)))./w;

xc = (V*(c_1.*(cos(w*(1:A)*dt)) + c_2.*(sin(w*(1:A)*dt))))./sm;
sc = (V*(c_1.*w.*(-sin(w*(1:A)*dt)) + c_2.*w.*(cos(w*(1:A)*dt))))./sm;

subplot(2,3,1);
plot(xc'+repelem(5*(1:3),A,1),'LineWidth',2);
legend("x_1 - time","x_2 - time","x_3 - time");
xlabel('time');
ylabel('position');
grid on;
title("Eigenvalues/vectors method - 'position-time'");

subplot(2,3,4);
plot(sc'+repelem(5*(1:3),A,1),'LineWidth',2);
legend("v_1 - time","v_2 - time","v_3 - time");
xlabel('time');
ylabel('velocity');
grid on;
title("Eigenvalues/vectors method - 'velocity-time'");


%Solution with Verlet algorithm
M = [
-k1/m1 k1/m1 0
k1/m2 -(k1+k2)/m2 k2/m2
0 k2/m3 -k2/m3
];

for t = 1 : A - 1
    a(:,t) = M*x(:,t);
    x(:,t+1) = x(:,t) + v(:,t)*dt + a(:,t)*dt^2/2;
    a(:,t+1) = M*x(:,t+1);
    v(:,t+1) = v(:,t) + (a(:,t) + a(:,t+1))/2*dt;
end

subplot(2,3,2)
plot(x'+repelem(5*(1:3),A,1),'LineWidth',2);
legend("x_1 - time","x_2 - time","x_3 - time");
xlabel('time');
ylabel('position');
grid on;
title("Verlet algorithm method 'position-time'");

subplot(2,3,5)
plot(v'+repelem(5*(1:3),A,1),'LineWidth',2);
legend("v_1 - time","v_2 - time","v_3 - time");
xlabel('time');
ylabel('velocity');
grid on;
title("Verlet algorithm method 'velocity-time'");

%Comparison
subplot(2,3,3)
plot(xc'+repelem(5*(1:3),A,1),"--",'LineWidth',2);
hold on;
plot(x'+repelem(5*(1:3),A,1),":",'LineWidth',2);
title("Comparison for 'position-time'")
xlabel('Time');
ylabel('Position');
grid on;

subplot(2,3,6)
plot(sc'+repelem(5*(1:3),A,1),"--",'LineWidth',2);
hold on;
plot(v'+repelem(5*(1:3),A,1),":",'LineWidth',2);
title("Comparison for 'velocity-time'")
xlabel('Time');
ylabel('Velocity');
grid on;