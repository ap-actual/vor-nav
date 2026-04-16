% This script will generate truth IMU data, then apply accelerometer spec
% errors, then apply gyro spec errors
clear all
close all
clc

load("filtered_flights.mat");
flightNames = fieldnames(flightData);
badFlights = {'icao_abaadd_F1';'icao_abaadd_F2';'icao_abd230_F1';...
    'icao_aaa28d_F1';'icao_aaad6d_F1';'icao_abef2f_F1';'icao_abdce1_F1';...
    'icao_abe0df_F1'};
flightNames(ismember(flightNames,badFlights)) = [];
imuTruth = struct();

for k = 1:numel(flightNames)
    thisName = flightNames{k};
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
    imuTruth.(thisName) = buildImuFromFlightData( ...
        lat, lon, velocity, heading, vertrate, geoalt, time);

    % build imuSpec
    imuSpec.(thisName) = buildImuSpec(); 
    
    % create imuMeasurements
    imuMeasured.(thisName) = addImuErrors(imuTruth.(thisName),imuSpec.(thisName));

end

%% RANDOM PLOTS
% %% MAKE PLOTS TO TEST IMU DATA INTEGRATION
% % test integrating to see if get same original NED velocity
% velNED = integrateSpecificForceNav(imuTruth.(thisName));
% vMag = vecnorm(velNED, 2, 2);
%
% figure()
% plot(imuTruth.(thisName).time,vMag)
% hold on; grid on
% plot(time,velocity* 0.514444)
% legend("Integrated Vel", "Original Vel")
% title("Integrate Specific Force vs. Velocity from Flight " + thisName)
% 
% figure()
% % make plot of truth and measured
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyTruth(:,1))
% hold on; grid on
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyMeas(:,1));
% legend("Truth", "Measured")
% title("X axis specF")
% 
% figure()
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyTruth(:,2))
% hold on; grid on
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyMeas(:,2));
% legend("Truth", "Measured")
% title("Y axis specF")
% 
% figure()
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyTruth(:,3))
% hold on; grid on
% plot(imuMeasured.(thisName).time, imuMeasured.(thisName).specificForceBodyMeas(:,3));
% legend("Truth", "Measured")
% title("Z axis specF")


