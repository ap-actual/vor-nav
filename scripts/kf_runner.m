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
dt = trajImuData.time(2) - trajImuData.time(1);

% Initial Conditions
% x0
lla0 = [trajImuData.lat(1), trajImuData.lon(1), trajImuData.alt(1)];
initialUTC = datetime(2000, 1, 1, 0, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'dd-MMM-uuuu HH:mm:ss.SSS');
pos0 = lla2eci(lla0, datevec(initialUTC))';
velNED0 = trajImuData.velocityNED(1,:);
[velECEFX0, velECEFY0, velECEFZ0] = ned2ecefv(velNED0(1), velNED0(2), velNED0(3), lla0(1), lla0(2));
velECEF0 = [velECEFX0, velECEFY0, velECEFZ0];
[~, vel0] = ecef2eci(initialUTC, zeros(1,3), velECEF0);
orn0 = trajImuData.eulerAngles(1,:)';
x0 = [pos0; vel0; orn0;... % dynamic states, X0 will be set depending on trajectory
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

   currentUTC = initialUTC + seconds(dt * i);
   sf = trajImuData.specificForceBodyMeas(i,:);
   bodyRates = trajImuData.gyroRatesMeas(i,:);

    % -- KF PROPAGATION --
    x(:,i+1) = dynamics(x(:,i), sf, bodyRates, dt, currentUTC);
    F = jacobian_f(x(:,i), sf, bodyRates, dt, currentUTC);
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
function xnext = dynamics(x,a,td,dt,UTC)

    tB2ECI = calculateBody2ECI(x, UTC);

    dx = zeros(21,1);

    % Velocity into position
    dx(1:3) = x(4:6);

    % Specific force into velocity; subtracting out our expected bias and
    % scale factor
    dx(4:6) = tB2ECI * (a' - x(16:18) .* 9.80665e-3)./(1 + x(10:12) .* 9.80665e-3);

    % Angular rates into orientation; subtracting out our expected bias and
    % scale factor
    dx(7:9) = (td' - x(19:21))./(1 + x(13:15) .* 1e-3);

    % Don't propogate bias and scale factor
    dx(10:21) = zeros(12,1);

    % propogate
    xnext = x + dt*dx;

end

% Jacobian
function F = jacobian_f(x,a,td,dt,UTC)

    n = length(x);
    eps = 1e-6;
    F = zeros(n);
    
    % propogate to next state
    f0 = dynamics(x,a,td,dt,UTC);
    
    % perturbations
    for i=1:n
        xp = x;
        xp(i) = xp(i) + eps;
        
        fp = dynamics(xp,a,td,dt,UTC);
        
        F(:,i) = (fp-f0)/eps;
    end
end

% tBody2ECI
function tB2ECI = calculateBody2ECI(x, UTC)
    % body to NED
    phi = x(7); sphi = sind(phi); cphi = cosd(phi);
    theta = x(8); stheta = sind(theta); ctheta = cosd(theta);
    psi = x(9); spsi = sind(psi); cpsi = cosd(psi);
    
    tB2NED = [ctheta * cpsi, ctheta * spsi, -stheta;
        sphi * stheta * cpsi - cphi * spsi, sphi * stheta * spsi + cphi * cpsi, sphi * ctheta;
        cphi * stheta * cpsi + sphi * spsi, cphi * stheta * spsi - sphi * cpsi, cphi * ctheta];

    % NED to ECEF
    lla = eci2lla(x(1:3)', datevec(UTC));
    lat = lla(1); slat = sind(lat); clat = cosd(lat);
    lon = lla(2); slon = sind(lon); clon = cosd(lon);

    tNED2ECEF = [-slat * clon, -slat * slon, clat;
        -slon, clon, 0;
        -clat * clon, -clat * slon, -slat];

    % ECEF to ECI
    omegaE = 7.2921159e-5; % rad/s
    theta = omegaE * seconds(timeofday(UTC));

    tECEF2ECI = [...
         cos(theta), -sin(theta), 0;
         sin(theta),  cos(theta), 0;
         0,           0,          1];

    tB2ECI = tECEF2ECI * tNED2ECEF * tB2NED;

    
end

% Measurement
function zhat = h(x, navAids, vorMeasData)

    pos = x(1:3); % aircraft position

    n = numel(vorMeasData);
    zhat = zeros(n,1);

    for k = 1:n
        stationPos = navAids(vorMeasData(k).id).position;

        dx = pos(1) - stationPos(1);
        dy = pos(2) - stationPos(2);

        zhat(k) = atan2(dy, dx);
    end
end

% Measurement jacobian
function H = jacobian_h(x, navAids, vorMeasData)

    pos = x(1:3);
    n = numel(vorMeasData);

    H = zeros(n, 21);

    for k = 1:n
        stationPos = navAids(vorMeasData(k).id).position;

        dx = pos(1) - stationPos(1);
        dy = pos(2) - stationPos(2);

        r2 = dx^2 + dy^2;

        H(k,1) = -dy / r2;
        H(k,2) =  dx / r2;
    end
end