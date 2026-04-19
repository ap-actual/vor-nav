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
[imuMeasured, imuTruth] = createImuMeasurements(flightData, 'icao_abe79d_F1');

% load navaid data
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% VOR measurement rate
vorMeasRate = 1; % hz
vorMeasLast = 0;

nTimes = numel(imuMeasured.time);

% Initial Conditions

x = zeros(21,nTimes);
P = zeros(21,21,nTimes);

x(:,1) = x0;
P(:,:,1) = P0;

% run for-loop
for i = 1:numel(imuMeasured.time) 

    % -- KF PROPAGATION STEP HERE --
    x(:,i+1) = dynamics(x(:,i), a, td,dt,i);
    F = jacobian_f(x(:,i),a,td,dt,i);
    Phi = expm(F * dt);
    P(:,:,i+1) = Phi * P(:,:,i) * Phi' + Q;
    
    % wowee look at me I'm KF propagation math



    % if VOR meas happens, run VOR meas
    if imuMeasured.time(i) - vorMeasLast > 1/vorMeasRate
        
        % generate VOR measurement
        truthLla = [imuMeasured.lat(i), imuMeasured.lon(i), imuMeasured.alt(i)];

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
        vorMeasLast = imuMeasured.time(i);
    end

end


%% helper functions

function [imuMeasured, imuTruth] = createImuMeasurements(flightData, flightId)

    thisName = flightId;

    T = flightData.(thisName);   

    % pull required columns to build imu data
    time = T.time;
    lat = T.lat;
    lon = T.lon;
    velocity = T.velocity;
    heading = T.heading;
    vertrate = T.vertrate;
    onground = T.onground;
    geoalt = T.geoaltitude;

    
   % Prune the data from the flightData set to be only valid stuff (not
   % NANS, not giant gaps of time)
    validIdx = ~(isnan(time) | isnan(lat) | isnan(lon) | isnan(velocity) | ...
        isnan(heading) | isnan(vertrate) | isnan(geoalt));

    time = time(validIdx);
    lat = lat(validIdx);
    lon = lon(validIdx);
    velocity = velocity(validIdx);
    heading = heading(validIdx);
    vertrate = vertrate(validIdx);
    onground = onground(validIdx);
    geoalt = geoalt(validIdx);

    % look for large gaps of time 
    time = time - time(1);
    timeDiff = diff(time);
    tol = 50;  
    flightBounds = find(timeDiff > tol);

    startIdx = flightBounds(1) + 1;
    endIdx = length(time);
    % keep only valid segments
    time = time(startIdx:endIdx);
    lat = lat(startIdx:endIdx);
    lon = lon(startIdx:endIdx);
    velocity = velocity(startIdx:endIdx);
    heading = heading(startIdx:endIdx);
    vertrate = vertrate(startIdx:endIdx);
    onground = onground(startIdx:endIdx);
    geoalt = geoalt(startIdx:endIdx);

    % build truth IMU for this flight
    imuTruth.(thisName) = VorNav.buildImuFromFlightData( ...
        lat, lon, velocity, heading, vertrate, geoalt, time);

    % build imuSpec
    imuSpec.(thisName) = VorNav.buildImuSpec(); 
    
    % create imuMeasurements
    imuMeasured = VorNav.addImuErrors(imuTruth.(thisName),imuSpec.(thisName));

    % normalize IMU measurement times
    imuMeasured.time = imuMeasured.time - imuMeasured.time(1);

end

function xnext = dynamics(x,a,td,dt,i)
    dx = zeros(21,1);
    dx(1:3) = x(4:6);
    dx(4:6) = (a(:,i) - x(16:18).*9.80665e-3)./(1+ x(10:12).*9.80665e-3);
    dx(7:9) = (td(:,i)-x(19:21))./(1+x(13:15).*1e-3);
    dx(10:21) = zeros(12,1);
    
    xnext = x + dt*dx;

end