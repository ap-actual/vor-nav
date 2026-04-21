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

R = VOR_1_SIGMA^2;

%% generate truth IMU and VOR measurement data

% clear everything except flight data (in case previous step has been skipped)

% create IMU measurement 
trajImuData= VorNav.createImuMeasurements(flightData, 'icao_abe79d_F1');

% load navaid data
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% find vors along flight path
visibleVorIdents = VorNav.findVorsAlongFlightpath(trajImuData, navAids);

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
pos0ECEF = lla2ecef(lla0);
velNED0 = trajImuData.velocityNED(1,:) .* 0.514444;
[velECEFX0, velECEFY0, velECEFZ0] = ned2ecefv(velNED0(1), velNED0(2), velNED0(3), lla0(1), lla0(2));
velECEF0 = [velECEFX0, velECEFY0, velECEFZ0];
[~, vel0] = ecef2eci(initialUTC, pos0ECEF, velECEF0);
orn0 = rad2deg(trajImuData.eulerAngles(1,:)');
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
   tB2ECI = calculateBody2ECI(x(:,i), currentUTC);
   sf = trajImuData.specificForceBodyTruth(i,:);
   bodyRates = trajImuData.gyroRatesTruth(i,:);

    % -- KF PROPAGATION --
    x(:,i+1) = dynamics(x(:,i), sf, bodyRates, dt, tB2ECI);
    F = jacobian_f(x(:,i), sf, bodyRates, dt, tB2ECI);
    Phi = expm(F * dt);
    P(:,:,i+1) = Phi * P(:,:,i) * Phi' + Q;

    trueLLA = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];
    trueECI = lla2eci(trueLLA, datevec(currentUTC));

    % if VOR meas happens, run VOR meas
    if trajImuData.time(i) - vorMeasLast > 1/vorMeasRate
        
        % generate VOR measurement
        truthLla = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];

        % create vor meas
        vorMeasData = VorNav.vorMeas(truthLla, "lla", navAids);

        % generate VOR meas with noise
        vorNoise = VOR_1_SIGMA * randn(1,numel([vorMeasData.bearing_deg]));
        vorMeasBearing = [vorMeasData.bearing_deg];
        vorMeasBearing = vorNoise + vorMeasBearing;

        % measurement update
        % plot setup for debug
        % llaStateEst = eci2lla(x(1:3,i+1)', datevec(currentUTC));
        % geoscatter(llaStateEst(1), llaStateEst(2));
        % 
        % z = [deg2rad(vorMeasBearing'); zeros(25,1)];
        % H = jacobian_h(x(:,i+1), vorMeasData, visibleVorIdents, currentUTC);
        % 

        H = predictedBearingJacobian(x(:,i+1), vorMeasData, visibleVorIdents, currentUTC);

        % compute HX
        hx = predictedBearings(trueECI', vorMeasData, visibleVorIdents, currentUTC);
        z_vor = deg2rad(vorMeasBearing');

        innovation = mod(z_vor - hx + pi, 2*pi) - pi;  % wraps to (-pi, pi]

        % check for 3 sigma innovation
        badInd = abs(innovation) >= 3 * deg2rad(sqrt(R));
        if all(badInd)
            fprintf("All measurements thrown out at t = %d seconds!\r\n", trajImuData.time(i))
            continue
        end
        if any(badInd)
            H(badInd,:) = zeros(sum(badInd), 21);
        end

        K = P(:,:,i+1) * H' / (H * P(:,:,i+1) * H' + R);

        x(:, i+1) = x(:, i+1) + K * innovation;

        P(:,:,i+1) = (eye(numel(21)) - K * H) *  P(:,:,i+1);
        
        fprintf("Successful measurement update at t = %d seconds!\r\n", trajImuData.time(i))

        % reset last meas time
        vorMeasLast = trajImuData.time(i);
    end

end


%% helper functions
% Dynamics Function
function xnext = dynamics(x,a,td,dt,tB2ECI)

    dx = zeros(21,1);

    % Velocity into position
    dx(1:3) = x(4:6);

    % Specific force into velocity; subtracting out our expected bias and
    % scale factor
    % dx(4:6) = tB2ECI * (a' - x(16:18) .* 9.80665e-3)./(1 + x(10:12) .* 9.80665e-3);
    dx(4:6) = tB2ECI * a';

    % Angular rates into orientation; subtracting out our expected bias and
    % scale factor
    % dx(7:9) = (td' - x(19:21))./(1 + x(13:15) .* 1e-3);
    dx(7:9) = rad2deg(td');

    % Don't propogate bias and scale factor
    dx(10:21) = zeros(12,1);

    % propogate
    xnext = x + dt*dx;

end

% Jacobian
function F = jacobian_f(x,a,td,dt,tB2ECI)

    n = length(x);
    eps = 1e-6;
    F = zeros(n);
    
    % propogate to next state
    f0 = dynamics(x,a,td,dt,tB2ECI);
    
    % perturbations
    for i=1:n
        xp = x;
        xp(i) = xp(i) + eps;
        
        fp = dynamics(xp,a,td,dt,tB2ECI);
        
        F(:,i) = (fp-f0)/eps;
    end
end

% tBody2ECI
function tB2ECI = calculateBody2ECI(x, UTC)
    % body to NED
    tB2NED = angle2dcm(deg2rad(x(7)),deg2rad(x(8)),deg2rad(x(9)),'ZYX');

    lla = eci2lla(x(1:3)', datevec(UTC));
    tECEF2NED = dcmecef2ned(lla(1), lla(2));

    tECI2ECEF = dcmeci2ecef('IAU-2000/2006', datevec(UTC));

    tB2ECI = tECI2ECEF' * tECEF2NED' * tB2NED;
    % body to NED
    % phi = x(7); sphi = sind(phi); cphi = cosd(phi);
    % theta = x(8); stheta = sind(theta); ctheta = cosd(theta);
    % psi = x(9); spsi = sind(psi); cpsi = cosd(psi);
    % 
    % tB2NED = [ctheta * cpsi, ctheta * spsi, -stheta;
    %     sphi * stheta * cpsi - cphi * spsi, sphi * stheta * spsi + cphi * cpsi, sphi * ctheta;
    %     cphi * stheta * cpsi + sphi * spsi, cphi * stheta * spsi - sphi * cpsi, cphi * ctheta];
    % 
    % % NED to ECEF
    % lla = eci2lla(x(1:3)', datevec(UTC));
    % lat = lla(1); slat = sind(lat); clat = cosd(lat);
    % lon = lla(2); slon = sind(lon); clon = cosd(lon);
    % 
    % tNED2ECEF = [-slat * clon, -slat * slon, clat;
    %     -slon, clon, 0;
    %     -clat * clon, -clat * slon, -slat];
    % 
    % % ECEF to ECI
    % omegaE = 7.2921159e-5; % rad/s
    % theta = omegaE * seconds(timeofday(UTC));
    % 
    % tECEF2ECI = [...
    %      cos(theta), -sin(theta), 0;
    %      sin(theta),  cos(theta), 0;
    %      0,           0,          1];
    % 
    % tB2ECI = tECEF2ECI * tNED2ECEF * tB2NED;

    
end

% Measurement jacobian
% function H = jacobian_h(x, vorMeasData, visibleVorIdents, UTC)
% 
%     pos = x(1:3);  % 3x1 ECI position
%     n = numel(visibleVorIdents);
% 
%     H = zeros(n, 21);
% 
%     visibleIdx = find(contains(string(vertcat(vorMeasData.ident)), visibleVorIdents));
% 
%     for k = 1:numel(visibleIdx)
% 
%         stationLLA = [vorMeasData(visibleIdx(k)).lat, vorMeasData(visibleIdx(k)).lon, 0];
%         stationPosECI = lla2eci(stationLLA, datevec(UTC));  % 1x3
% 
%         datevecUTC = datevec(UTC);
% 
%         % ECI to ECEF via MATLAB built-in (returns 3x3 rotation matrix)
%         tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevecUTC);
% 
%         % ECEF to NED via MATLAB built-in
%         lat = stationLLA(1);  % degrees
%         lon = stationLLA(2);  % degrees
%         tECEF2NED = dcmecef2ned(lat, lon);  % built-in from Aerospace Toolbox
% 
%         tECI2NED = tECEF2NED * tECI2ECEF_mat;
% 
%         % FIX #4: transpose stationPosECI to match 3x1 pos
%         dr = pos - stationPosECI';   % 3x1
% 
%         drNED = tECI2NED * dr;       % 3x1, no transpose needed now
% 
%         dx = drNED(1);   % North
%         dy = drNED(2);   % East
% 
%         r2 = dx^2 + dy^2;
% 
%         % Hvec: d(bearing)/d(posECI), bearing = atan2(East, North)
%         Hvec = [-dy/r2, dx/r2, 0] * tECI2NED;  % 1x3
%         H(k, 1:3) = Hvec;
%     end
% end

function hx = predictedBearings(x, vorMeasData, visibleVorIdents, UTC)
    pos = x(1:3);
    visibleIdx = find(contains(string(vertcat(vorMeasData.ident)), visibleVorIdents));
    n = numel(visibleIdx);
    hx = zeros(n, 1);

    tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevec(UTC));

    for k = 1:numel(visibleIdx)
        stationLLA = [vorMeasData(k).lat, vorMeasData(k).lon, 0];
        stationPosECI = lla2eci(stationLLA, datevec(UTC))';

        lat = stationLLA(1);
        lon = stationLLA(2);
        tECEF2NED = dcmecef2ned(lat, lon);
        tECI2NED = tECEF2NED * tECI2ECEF_mat;

        dr = pos - stationPosECI;
        drNED = tECI2NED * dr;

        % atan2 returns (-pi, pi]; wrap to [0, 2pi)
        hx(k) = mod(atan2(drNED(2), drNED(1)), 2*pi);
    end
end

function H = predictedBearingJacobian(x, vorMeasData, visibleVorIdents, UTC)

    pos = x(1:3);
    visibleIdx = find(contains(string(vertcat(vorMeasData.ident)), visibleVorIdents));
    n = numel(visibleIdx);

    H = zeros(n, numel(x)); % assumes x may be larger than 3

    tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevec(UTC));

    for k = 1:n

        stationLLA = [vorMeasData(k).lat, vorMeasData(k).lon, 0];
        stationPosECI = lla2eci(stationLLA, datevec(UTC))';

        lat = stationLLA(1);
        lon = stationLLA(2);
        tECEF2NED = dcmecef2ned(lat, lon);
        T = tECEF2NED * tECI2ECEF_mat;

        dr = pos - stationPosECI;
        drNED = T * dr;

        N = drNED(1);
        E = drNED(2);

        denom = (N^2 + E^2);

        % Jacobian wrt position
        dHdpos = (-E * T(1,:) + N * T(2,:)) / denom;

        H(k,1:3) = dHdpos;
    end
end
