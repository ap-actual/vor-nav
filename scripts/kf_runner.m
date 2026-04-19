%% Kalman Filter Runner
% this script is designed to be RUN FROM THE SCRIPTS FOLDER

%% clean up and setup path
% add analysis scripts to path and setup path
clear,clc
cd ..
setup_path
load("dataSets\filtered_flights.mat");
cd scripts\

%% generate truth IMU and VOR measurement data

% clear everything except flight data (in case previous step has been skipped)
clearvars -except flightData

[imuMeasured, imuTruth] = buildGoldenDataset(flightData, 'icao_abe79d_F1')






%% helper functions

function [imuMeasured, imuTruth] = buildGoldenDataset(flightData, flightId)

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
    imuMeasured.(thisName) = VorNav.addImuErrors(imuTruth.(thisName),imuSpec.(thisName));

end