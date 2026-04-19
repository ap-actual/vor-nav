%% EKF Implementation
% Assuming ISM330DHCX MEMS IMU
% +- 4 g Acceleration Mode
% +- 125 dps Angular Rate Mode
% State, Frame, Units
% Position X, Y, Z, ECI, m
% Velocity X, Y, Z, ECI, m/s
% Orientation U, V, W, ECI, deg
% Linear Acceleration Sensitivity X, Y, Z, Accelerometer, mg/LSB 
% Angular Rate Sensitivity X, Y, Z, Gyro, mdps/LSB
% Linear acceleration zero-g level offset accuracy, X, Y, Z, Accelerometer, mg
% Angular rate zero-rate level, X, Y, Z, Gyro, dps

% Process Noise
% Rn 5-8 mdps/sqrt(hz)
% ARW = 0.21-0.34 deg/sqrt(hr)
% An 60-100 ug/sqrt(hz)
a = zeros(3,20);
a(1,:) = ones(1,20) .* 9.81;
td = zeros(3,20);
dt = 1/100;

% x0
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0;... % dynamic states, X0 will be set depending on trajectory
    0.122; 0.122; 0.122; ...        % Accel Sensitivies, based on document
    4.375; 4.375; 4.375;...         % Gyro Sensitivies, based on document
    0; 0; 0; 0; 0; 0];                          % zero-g level offsets 

% P0
P0 = diag([10000; 10000; 10000; 100; 100; 100; 180; 180; 180;... % Prior free (?) dynamic states
    6.4000e-07; 6.4000e-07; 6.4000e-07;...
    8.5069e-04; 8.5069e-04; 8.5069e-04;...
    469.4444; 469.4444; 469.4444;...
    1; 1; 1]);

x = zeros(21,size(a,2)+1);
P = zeros(21,21,size(a,2)+1);
x(:,1) = x0;
P(:,:,1) = P0;

% Q 
Qc = diag([0; 0; 0; 3.4621e-07; 3.4621e-07; 3.4621e-07; 1.2250e-05; 1.2250e-05; 1.2250e-05; zeros(12,1)]);
Q = Qc * dt;

% H unknown

% R unknown
update = 0;

for i = 1:size(a,2)
    x(:,i+1) = dynamics(x(:,i), a, td,dt,i);
    F = jacobian_f(x(:,i),a,td,dt,i);
    Phi = expm(F * dt);
    P(:,:,i+1) = Phi * P(:,:,i) * Phi' + Q;
    if update
        H = jacobian_h(x, a, td, dt,i);

    end
end

%% =========================
%% Dynamics Function
%% =========================

function xnext = dynamics(x,a,td,dt,i)
dx = zeros(21,1);
dx(1:3) = x(4:6);
dx(4:6) = (a(:,i) - x(16:18).*9.80665e-3)./(1+ x(10:12).*9.80665e-3);
dx(7:9) = (td(:,i)-x(19:21))./(1+x(13:15).*1e-3);
dx(10:21) = zeros(12,1);

xnext = x + dt*dx;

end

%% =========================
%% Jacobian (numerical)
%% =========================

function F = jacobian_f(x,a,td,dt,k)

n=length(x);
eps=1e-6;
F=zeros(n);

f0=dynamics(x,a,td,dt,k);

for i=1:n
xp=x;
xp(i)=xp(i)+eps;

fp=dynamics(xp,a,td,dt,k);

F(:,i)=(fp-f0)/eps;
end

end

function H = jacobian_h(x,a,td,dt)

n=length(x);
eps=1e-6;
F=zeros(n);

f0=dynamics(x,a,td,dt);

for i=1:n
xp=x;
xp(i)=xp(i)+eps;

fp=dynamics(xp,a,td,dt);

F(:,i)=(fp-f0)/eps;
end

end