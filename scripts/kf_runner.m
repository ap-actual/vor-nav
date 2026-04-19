%% Kalman Filter Runner
% this script is designed to be RUN FROM THE SCRIPTS FOLDER

%% clean up and setup path
% add analysis scripts to path and setup path
clear,clc
cd ..
setup_path
load("dataSets\filtered_flights.mat");
cd scripts\

%% set up script constants
clearvars -except flightData

VOR_1_SIGMA = 2; % 2-deg 1-sigma


%% generate truth IMU and VOR measurement data

% clear everything except flight data (in case previous step has been skipped)

% create IMU measurement 
trajImuData= VorNav.createImuMeasurements(flightData, 'icao_abe79d_F1');

% load navaid data
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% VOR measurement rate
vorMeasRate = 1; % hz
vorMeasLast = 0;

nTimes = numel(trajImuData.time);

% Initial Conditions
% x0
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0;... % dynamic states, X0 will be set depending on trajectory
    0.122; 0.122; 0.122; ...        % Accel Sensitivies, based on document
    4.375; 4.375; 4.375;...         % Gyro Sensitivies, based on document
    0; 0; 0; ...                    % Zero G Accel
    0; 0; 0];                       % Zero G Gyro

% P0
P0 = diag([1e6; 1e6; 1e6; 100; 100; 100; 180; 180; 180;... % Dynamics prior free(?)
    6.4000e-07; 6.4000e-07; 6.4000e-07;...                 % Accel Sensitivites
    8.5069e-04; 8.5069e-04; 8.5069e-04;...                 % Gyro Sensitivites
    469.4444; 469.4444; 469.4444;...                       % Zero G Accel
    1; 1; 1]);                                             % Zero G Gyro

% Q 
Qc = diag([0; 0; 0; 3.4621e-07; 3.4621e-07; 3.4621e-07; 1.2250e-05; 1.2250e-05; 1.2250e-05; zeros(12,1)]);
Q = Qc * dt;

x = zeros(21,nTimes);
P = zeros(21,21,nTimes);

x(:,1) = x0;
P(:,:,1) = P0;

% run for-loop
for i = 1:numel(trajImuData.time) 

    % -- KF PROPAGATION --
    x(:,i+1) = dynamics(x(:,i), a, td,dt,i);
    F = jacobian_f(x(:,i),a,td,dt,i);
    Phi = expm(F * dt);
    P(:,:,i+1) = Phi * P(:,:,i) * Phi' + Q;

    % if VOR meas happens, run VOR meas
    if trajImuData.time(i) - vorMeasLast > 1/vorMeasRate
        
        % generate VOR measurement
        truthLla = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];

        % create vor meas
        vorMeasData = VorNav.vorMeas(truthLla, "lla", navAids);

        % -- KF MEASUREMENT UPDATE HERE --
        % note: bearing measurements from VorNav are truth meas, will need
        % to add random draw to them

        % generate VOR meas with noise
        vorNoise = VOR_1_SIGMA * randn(1,numel([vorMeasData.bearing_deg]));
        vorMeasBearing = [vorMeasData.bearing_deg];
        vorMeasBearing = vorNoise + vorMeasBearing;

        % UPDATE LOGIC GOES HERE



        % reset last meas time
        vorMeasLast = trajImuData.time(i);
    end

end


%% helper functions
% Dynamics Function
function xnext = dynamics(x,a,td,dt,i)
    dx = zeros(21,1);
    dx(1:3) = x(4:6);
    dx(4:6) = (a(:,i) - x(16:18).*9.80665e-3)./(1+ x(10:12).*9.80665e-3);
    dx(7:9) = (td(:,i)-x(19:21))./(1+x(13:15).*1e-3);
    dx(10:21) = zeros(12,1);
    
    xnext = x + dt*dx;

end

% Jacobian
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