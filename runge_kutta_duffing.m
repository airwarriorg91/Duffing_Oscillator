%solving the duffing equation using 4th Order Runge Kutta
m = 0.2;
mu = 5;
cv = 1;
f = 2048; % sampling frequency 


%derivative function f'(x) = Ax + Bgamma
%Initial conditions
x1_0 = 0.05; %amplitude
x2_0 = 0; %zero velocity at extreme position
Y_0 = [x1_0;x2_0];  

h = 1/f; %steps for Runge Kutta
t = 0:h:5; %time 
y = zeros(2, length(t)); % preallocating the output
y(:,1)=Y_0;
G = zeros(1, length(t));
k = zeros(1, length(t));
for i = 1:(length(t)-1)
    y_i = y(:,i);
    k(:,i)=(K(y_i(1))*y_i(1))+(m*9.8);
    A = [0 1; -K(y_i(1)) -(cv + (mu*y_i(2)*signum(y_i(2))))];
    B = [0 ; Gamma(m)];
    G(:,i)=Gamma(m);
    m1 = A*y_i + B;
    m2 = A*(y_i + (m1*0.5*h)) + B;
    m3 = A*(y_i + (m2*0.5*h)) + B;
    m4 = A*(y_i + (m3*h)) + B;
    y(:,i+1)= y_i + ((m1+(2*m2)+(2*m3)+m4)/6)*h;
end

f = figure;
f.Units='inches';

f.Position(3:4)=[6 11];
subplot(3,1,1)
plot(y(1,:),y(2,:))
grid on;
xlabel("Displacement (m)")
ylabel("Velocity (m/s)")
title("Displacement v/s Velocity Curve for a Double Well Duffing Oscillator")

subplot(3,1,2)
plot(t,y(1,:))
grid on;
ylabel("Displacement (m)")
xlabel("Time (s)")
title("Displacement Curve for a Double Well Duffing Oscillator")

subplot(3,1,3)
plot(t,y(2,:))
ylabel("Velocity (m/s)")
xlabel("Time (s)")
grid on;
title("Velocity Curve for a Double Well Duffing Oscillator")

g = figure();
plot(y(1,:),k)
grid on;
title("Restoring Force of the Double Well Duffing Oscillator system")
xlabel("Displacement (m)")
ylabel("K (N)")

G_prob = gauss_distribution(G,0,50);
h = figure();
h.Units='inches';
h.Position(3:4)=[8 5];
scatter(G,G_prob,1);
grid on;
title("Input Perturbation Force (Zero Mean Random Gaussian Function)")
ylabel("Probability")
xlabel("Force (N)")
