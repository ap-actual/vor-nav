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

VOR_1_SIGMA = deg2rad(2); % 2-deg 1-sigma

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
% gather info for x0
lla0 = [trajImuData.lat(1), trajImuData.lon(1), trajImuData.alt(1)];
initialUTC = datetime(2000, 1, 1, 0, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'dd-MMM-uuuu HH:mm:ss.SSS');
pos0 = lla2eci(lla0, datevec(initialUTC))';
pos0ECEF = lla2ecef(lla0);
velNED0 = trajImuData.velocityNED(1,:) .* 0.514444;
[velECEFX0, velECEFY0, velECEFZ0] = ned2ecefv(velNED0(1), velNED0(2), velNED0(3), lla0(1), lla0(2));
velECEF0 = [velECEFX0, velECEFY0, velECEFZ0];
[~, vel0] = ecef2eci(initialUTC, pos0ECEF, velECEF0);
orn0 = rad2deg(trajImuData.eulerAngles(1,:)');

% x0
x0 = [pos0; vel0; orn0];

% P0
P0 = diag([1e8; 1e8; 1e8; 10; 10; 10; pi/4; pi/4; pi/4]);

% Q
% Qc = diag([0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1])*10000;
% Q = Qc * dt;
Qc = diag([100; 100; 100; 1; 1; 1; 0.01; 0.01; 0.01])*1000;
Q = Qc * dt;

% pre-allocate state
x = zeros(9,nTimes);
P = zeros(9,9,nTimes);

% set initial state
x(:,1) = x0;
P(:,:,1) = P0;

% run for-loop
for i = 1:numel(trajImuData.time) 

   currentUTC(i) = initialUTC + seconds(dt * i);
   tB2ECI = calculateBody2ECI(x(:,i), currentUTC(i));
   sf = trajImuData.specificForceBodyTruth(i,:);
   bodyRates = trajImuData.gyroRatesTruth(i,:);

    % -- KF PROPAGATION --
    x(:,i+1) = dynamics(x(:,i), sf, bodyRates, dt, tB2ECI);
    F = jacobian_f(x(:,i), sf, bodyRates, dt, tB2ECI);
    Phi = F;
    P(:,:,i+1) = Phi * P(:,:,i) * Phi' + Q;

    % after computing F/Phi, before propagating P
    eigsPhi = eig(Phi);

    trueLLA = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];
    trueECI = lla2eci(trueLLA, datevec(currentUTC(i)));

    % if VOR meas happens, run VOR meas
    if trajImuData.time(i) - vorMeasLast > 1/vorMeasRate
        
        % generate VOR measurement
        truthLla = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];

        % create vor meas
        vorMeasData = VorNav.vorMeas(truthLla, "lla", navAids);
        nVor = numel(vorMeasData);
        R = VOR_1_SIGMA^2 * eye(nVor);

        % generate VOR meas with noise
        vorNoise = VOR_1_SIGMA * randn(1,numel([vorMeasData.bearing_deg]));
        vorMeasBearing = [vorMeasData.bearing_deg];
        vorMeasBearing = vorNoise + vorMeasBearing;

        % compute jacobian
        H = predictedBearingJacobian(x(:,i+1), vorMeasData, visibleVorIdents, currentUTC(i));

        % compute HX
        hx = predictedBearings(x(1:3,i+1), vorMeasData, visibleVorIdents, currentUTC(i));
        z_vor = deg2rad(vorMeasBearing');

        innovation = mod(z_vor - hx + pi, 2*pi) - pi;  % wraps to (-pi, pi]
        % disp(innovation)

        S = H * P(:,:,i+1) * H' + R;
        S_diag = diag(S);  % individual innovation variances

        % find which measurements pass the gate
        valid = false(nVor, 1);
        for j = 1:nVor
            norm_innov_sq = innovation(j)^2 / S_diag(j);  % scalar chi2, 1-DOF
            if norm_innov_sq < chi2inv(0.9973, 1)  % = 8.99
                valid(j) = true;
            else
                fprintf("Bearing %.0f gated at t = %.1f s (%.2f sigma)\n", ...
                    j, trajImuData.time(i), sqrt(norm_innov_sq));
            end
        end

        % only update if at least one measurement passed
        if any(valid)
            H_v = H(valid, :);
            R_v = R(valid, valid);
            innovation_v = innovation(valid);

            S_v = H_v * P(:,:,i+1) * H_v' + R_v;
            K = P(:,:,i+1) * H_v' / S_v;
            x(:,i+1) = x(:,i+1) + K * innovation_v;
            IKH = eye(9) - K * H_v;
            P(:,:,i+1) = IKH * P(:,:,i+1) * IKH' + K * R_v * K';
            P(:,:,i+1) = (P(:,:,i+1) + P(:,:,i+1)') / 2;
            fprintf("Measurement update at t = %f s (%d/%d bearings used)\n", ...
                trajImuData.time(i), sum(valid), nVor);

            eigvals = eig(P(:,:,i+1));
            fprintf("min eigenvalue of P: %.4e\n", min(eigvals));
        end
        
        % reset last meas time
        vorMeasLast = trajImuData.time(i);

        % ------ GEOPLOT UPDATE -------------------------------------------
        figure(1); clf;

        % plot all known VORs (dim)
        % allLats = [navAids.high(:).y, navAids.low(:).y];
        % allLons = [navAids.high(:).x, navAids.low(:).x];
        % geoplot(allLats, allLons, 'b.', 'MarkerSize', 6, 'DisplayName', 'All VORs');
        % hold on;

        % plot visible/active VORs (bright, labeled)
        if any(valid)
            activeLats = [vorMeasData(valid).lat];
            activeLons = [vorMeasData(valid).lon];
            geoplot(activeLats, activeLons, 'b^', 'MarkerSize', 10, ...
                'MarkerFaceColor', 'b', 'DisplayName', 'Active VORs');
            hold on;

            for j = 1:numel(activeLats)
                text(activeLats(j), activeLons(j), ...
                    sprintf('  %s', vorMeasData(j).ident), ...  % adjust field name to match your struct
                    'FontSize', 8, 'Color', 'blue');
            end
        end

        % convert current ECI state estimate to LLA for plotting
        estLLA = eci2lla(x(1:3, i+1)', datevec(currentUTC(i)));
        geoplot(estLLA(1), estLLA(2), 'g*', 'MarkerSize', 14, ...
            'MarkerFaceColor', 'g', 'DisplayName', 'KF Estimate');

        % plot truth for reference
        geoplot(trueLLA(1), trueLLA(2), 'r+', 'MarkerSize', 12, ...
            'LineWidth', 2, 'DisplayName', 'Truth');

        % -- VOR BEARING CONES --
        if any(valid)
            visIdx = find(valid);
            CONE_HALF_ANGLE_DEG = 2;

            for j = 1:numel(visIdx)
                vorLat = vorMeasData(visIdx(j)).lat;
                vorLon = vorMeasData(visIdx(j)).lon;

                % bearing from VOR to current position ESTIMATE (not raw measurement)
                dLat = estLLA(1) - vorLat;
                dLon = estLLA(2) - vorLon;
                bearingToEst = atan2(dLon, dLat);  % N=0, E=pi/2 convention

                % cone length = exact distance from VOR to estimate
                CONE_LENGTH_DEG = sqrt(dLat^2 + dLon^2);

                % left and right cone edges centered on bearing to estimate
                leftBearing  = bearingToEst - deg2rad(CONE_HALF_ANGLE_DEG);
                rightBearing = bearingToEst + deg2rad(CONE_HALF_ANGLE_DEG);

                % arc endpoint at aircraft distance
                arcAngles = linspace(leftBearing, rightBearing, 30);
                arcLats = vorLat + CONE_LENGTH_DEG * cos(arcAngles);
                arcLons = vorLon + CONE_LENGTH_DEG * sin(arcAngles);

                coneLats = [vorLat, arcLats, vorLat];
                coneLons = [vorLon, arcLons, vorLon];

                % outline (all MATLAB versions)
                geoplot(coneLats, coneLons, 'c-', 'LineWidth', 1, 'DisplayName', '');

                % % filled cone (R2023a+ only, comment out if unavailable)
                % geopatch(coneLats, coneLons, 'c', 'FaceAlpha', 0.12, ...
                %     'EdgeColor', 'cyan', 'LineWidth', 1, 'DisplayName', '');
            end

            % single legend entry
            geoplot(NaN, NaN, 'c-', 'LineWidth', 1.5, 'DisplayName', '2° VOR cone');
        end


        % -- covariance ellipse in LLA --
        % extract position covariance and project to lat/lon degrees
        P_pos = P(1:3, 1:3, i+1);          % 3x3 ECI position covariance (meters²)
        METERS_PER_DEG = 111139;

        % rough ECI->lat/lon covariance projection (diagonal approximation)
        P_latlon = P_pos(1:2, 1:2) / METERS_PER_DEG^2;  % 2x2 in degrees²

        % eigendecomposition for ellipse axes
        [V, D] = eig(P_latlon);
        angles = linspace(0, 2*pi, 100);
        nSigma = 3;
        unit_circle = [cos(angles); sin(angles)];
        ellipse_pts = V * (nSigma * sqrt(D)) * unit_circle;  % 2xN in degrees
        ellipse_lats = estLLA(1) + ellipse_pts(1, :);
        ellipse_lons = estLLA(2) + ellipse_pts(2, :);
        geoplot(ellipse_lats, ellipse_lons, 'y-', 'LineWidth', 1.5, ...
            'DisplayName', '3\sigma bounds');

        % formatting
        geobasemap('streets');   % or 'satellite', 'topographic', 'darkwater'
        legend('Location', 'best');
        title(sprintf('KF State  |  t = %.1f s  |  %d/%d bearings used', ...
            trajImuData.time(i), sum(valid), nVor));
        drawnow;
    end

end


%% helper functions
% Dynamics Function
function xnext = dynamics(x,a,bodyRates,dt,tB2ECI)
    % STATE VECTOR:
    % x(1:3) - ECI position (m)
    % x(4:6) - ECI velocity (m/s)
    % x(7:9) - Euler angles [phi, theta, psi] (deg) ZYX convention

    dx = zeros(9,1);

    % --- POSITION: velocity into position ---
    dx(1:3) = x(4:6);

    % --- VELOCITY: specific force rotated to ECI + gravity ---
    mu = 3.986004418e14;        % m^3/s^2
    r_vec = x(1:3);
    r_norm = norm(r_vec);
    g_ECI = -mu / r_norm^3 * r_vec;

    dx(4:6) = tB2ECI * a' + g_ECI;

    % --- ORIENTATION: Euler kinematic equations ---
    phi   = deg2rad(x(7));
    theta = deg2rad(x(8));

    p = bodyRates(1);
    q = bodyRates(2);
    r = bodyRates(3);

    % guard against gimbal lock (theta near +/- 90 deg)
    cos_theta = cos(theta);
    if abs(cos_theta) < 1e-6
        cos_theta = sign(cos_theta) * 1e-6;
    end

    phi_dot   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
    theta_dot =     q*cos(phi)            - r*sin(phi);
    psi_dot   =    (q*sin(phi)            + r*cos(phi)) / cos_theta;

    dx(7:9) = rad2deg([phi_dot; theta_dot; psi_dot]);

    % --- PROPAGATE ---
    xnext = x + dt * dx;


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
    n = numel(vorMeasData);  % iterate directly, no visibleIdx filtering needed
    hx = zeros(n, 1);
    tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevec(UTC));

    for k = 1:n
        stationLLA = [vorMeasData(k).lat, vorMeasData(k).lon, 0];
        stationPosECI = lla2eci(stationLLA, datevec(UTC))';
        lat = stationLLA(1);
        lon = stationLLA(2);
        tECEF2NED = dcmecef2ned(lat, lon);
        tECI2NED = tECEF2NED * tECI2ECEF_mat;
        dr = pos - stationPosECI;
        drNED = tECI2NED * dr;
        hx(k) = mod(atan2(drNED(2), drNED(1)), 2*pi);
    end
end


function H = predictedBearingJacobian(x, vorMeasData, visibleVorIdents, UTC)
    pos = x(1:3);
    n = numel(vorMeasData);  % iterate directly, no visibleIdx filtering needed
    H = zeros(n, numel(x));
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
        denom = max(N^2 + E^2, 1e-6);
        dHdpos = (-E * T(1,:) + N * T(2,:)) / denom;
        H(k, 1:3) = dHdpos;
    end
end
