%% Kalman Filter Runner - Integrated Simulation & Visualization
% This script runs the KF, saves results to a .mat file, and then animates.
clear; clc; close all;

%% 1. Setup and Data Loading
% cd .. % Uncomment if running from scripts folder
% setup_path
load("dataSets\filtered_flights.mat");
% cd scripts\

%% 2. Constants & Settings
VOR_1_SIGMA = deg2rad(2); 
R_noise      = VOR_1_SIGMA^2;
GATE_LIMIT   = 9.0; % 3-sigma gate

%% 3. Generate Truth & Measurements
trajImuData = VorNav.createImuMeasurements(flightData, 'icao_abe79d_F1');
navAids     = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

vorMeasRate = 1; % Hz
vorMeasLast = -inf;
nTimes      = numel(trajImuData.time);
dt          = trajImuData.time(2) - trajImuData.time(1);

%% 4. Initial Conditions
lla0       = [trajImuData.lat(1), trajImuData.lon(1), trajImuData.alt(1)];
initialUTC = datetime(2000, 1, 1, 0, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'dd-MMM-uuuu HH:mm:ss.SSS');
pos0ECEF   = lla2ecef(lla0);
pos0       = lla2eci(lla0, datevec(initialUTC))';
velNED0    = trajImuData.velocityNED(1,:) .* 0.514444; % knots to m/s
[vE_x, vE_y, vE_z] = ned2ecefv(velNED0(1), velNED0(2), velNED0(3), lla0(1), lla0(2));
[~, vel0]  = ecef2eci(initialUTC, pos0ECEF, [vE_x, vE_y, vE_z]);
orn0       = rad2deg(trajImuData.eulerAngles(1,:)');

% State: [Pos_ECI(3), Vel_ECI(3), Euler(3)]
x0 = [pos0; vel0; orn0];
P0 = diag([1e8; 1e8; 1e8; 10; 10; 10; pi/4; pi/4; pi/4]);
Qc = diag([100; 100; 100; 10; 10; 10; 0.01; 0.01; 0.01])*10;
Q  = Qc * dt;

%% 5. Pre-allocate Results Structure
results.time     = trajImuData.time;
results.x        = zeros(9, nTimes);
results.P        = zeros(9, 9, nTimes);
results.nUsed    = zeros(1, nTimes);
results.nVisible = zeros(1, nTimes);
results.utc      = datetime(zeros(nTimes, 6), 'TimeZone', 'UTC');

results.x(:,1)   = x0;
results.P(:,:,1) = P0;

%% 6. Main Filter Loop (Computation Only)
fprintf('--- Starting Kalman Filter ---\n');
tic;
for i = 1:nTimes-1
    results.utc(i) = initialUTC + seconds(dt * (i-1));
    
    % --- KF PROPAGATION ---
    tNED2ECI  = calculateNED2ECI(results.x(:,i), results.utc(i));
    sf        = trajImuData.specificForceNav(i,:);
    bodyRates = trajImuData.gyroRatesTruth(i,:);
    
    results.x(:,i+1)    = dynamics(results.x(:,i), sf, bodyRates, dt, tNED2ECI);
    F                   = jacobian_f(results.x(:,i), sf, bodyRates, dt, tNED2ECI);
    results.P(:,:,i+1)  = F * results.P(:,:,i) * F' + Q;
    
    % --- VOR MEASUREMENT UPDATE ---
    if trajImuData.time(i) - vorMeasLast >= 1/vorMeasRate
        truthLla    = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];
        vorMeasData = VorNav.vorMeas(truthLla, "lla", navAids);
        nVor        = numel(vorMeasData);
        results.nVisible(i) = nVor;
        
        vorNoise       = VOR_1_SIGMA * randn(1, nVor);
        vorMeasBearing = deg2rad([vorMeasData.bearing_deg] + vorNoise);
        
        nUsed = 0;
        for j = 1:nVor
            H_all  = predictedBearingJacobian(results.x(:,i+1), vorMeasData, results.utc(i));
            hx_all = predictedBearings(results.x(1:3,i+1), vorMeasData, results.utc(i));
            
            H_j          = H_all(j, :);
            hx_j         = hx_all(j);
            innovation_j = mod(vorMeasBearing(j) - hx_j + pi, 2*pi) - pi;
            S_j          = H_j * results.P(:,:,i+1) * H_j' + R_noise;
            
            % Gating
            if (innovation_j^2 / S_j) <= GATE_LIMIT
                K = results.P(:,:,i+1) * H_j' / S_j;
                results.x(:,i+1)   = results.x(:,i+1) + K * innovation_j;
                IKH                = eye(9) - K * H_j;
                results.P(:,:,i+1) = IKH * results.P(:,:,i+1) * IKH' + K * R_noise * K';
                nUsed = nUsed + 1;
            end
        end

        % --- TRUTH ALTITUDE INJECTION ---
        % Convert current ECI estimate to LLA, replace alt with truth, convert back
        est_lla = eci2lla(results.x(1:3, i+1)', datevec(results.utc(i)));
        truth_alt_m = trajImuData.alt(i+1) * 0.3048; % ft -> m (if alt is in feet)
        corrected_lla = [est_lla(1), est_lla(2), truth_alt_m];
        results.x(1:3, i+1) = lla2eci(corrected_lla, datevec(results.utc(i)))';

        results.nUsed(i) = nUsed;
        vorMeasLast = trajImuData.time(i);
        fprintf("t = %.2f\n", trajImuData.time(i));
    end
end
results.utc(end) = initialUTC + seconds(dt * (nTimes-1));
toc;

%% 7. Save Data
save('kf_results_icao_abe79d_F1_v2.mat', 'results', 'trajImuData', 'navAids');
fprintf('Data saved to kf_results.mat\n');


%% 8. Post-Processing Visualization (High-Speed Handle Graphics)
VorNav.PlotUtilies.plotStateEstimate(results.utc, results.x, results.P, trajImuData)

%% ========================================================================
%  Helper Functions
%  ========================================================================

function xnext = dynamics(x, a, bodyRates, dt, tNED2ECI)
    dx = zeros(9,1);
    dx(1:3) = x(4:6); % Velocity
    mu = 3.986004418e14;
    r_vec = x(1:3);
    g_ECI = -mu / norm(r_vec)^3 * r_vec;
    dx(4:6) = tNED2ECI * a' + g_ECI; % Accel
    
    phi = deg2rad(x(7)); theta = deg2rad(x(8));
    p = bodyRates(1); q = bodyRates(2); r = bodyRates(3);
    
    cT = max(abs(cos(theta)), 1e-6) * sign(cos(theta));
    phi_dot   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
    theta_dot = q*cos(phi) - r*sin(phi);
    psi_dot   = (q*sin(phi) + r*cos(phi)) / cT;
    dx(7:9) = rad2deg([phi_dot; theta_dot; psi_dot]);
    xnext = x + dt * dx;
end

function F = jacobian_f(x, a, bodyRates, dt, tNED2ECI)
    n = length(x); eps_val = 1e-6; F = zeros(n);
    f0 = dynamics(x, a, bodyRates, dt, tNED2ECI);
    for i = 1:n
        xp = x; xp(i) = xp(i) + eps_val;
        fp = dynamics(xp, a, bodyRates, dt, tNED2ECI);
        F(:,i) = (fp - f0) / eps_val;
    end
end

function tNED2ECI = calculateNED2ECI(x, UTC)
    lla = eci2lla(x(1:3)', datevec(UTC));
    tECEF2NED = dcmecef2ned(lla(1), lla(2));
    tECI2ECEF = dcmeci2ecef('IAU-2000/2006', datevec(UTC));
    tNED2ECI  = tECI2ECEF' * tECEF2NED';
end

function hx = predictedBearings(x, vorMeasData, UTC)
    pos = x(1:3); n = numel(vorMeasData); hx = zeros(n, 1);
    tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevec(UTC));
    for k = 1:n
        sLLA = [vorMeasData(k).lat, vorMeasData(k).lon, 0];
        sPosECI = lla2eci(sLLA, datevec(UTC))';
        T = dcmecef2ned(sLLA(1), sLLA(2)) * tECI2ECEF_mat;
        drNED = T * (pos - sPosECI);
        hx(k) = mod(atan2(drNED(2), drNED(1)), 2*pi);
    end
end

function H = predictedBearingJacobian(x, vorMeasData, UTC)
    pos = x(1:3); n = numel(vorMeasData); H = zeros(n, 9);
    tECI2ECEF_mat = dcmeci2ecef('IAU-2000/2006', datevec(UTC));
    for k = 1:n
        sLLA = [vorMeasData(k).lat, vorMeasData(k).lon, 0];
        sPosECI = lla2eci(sLLA, datevec(UTC))';
        T = dcmecef2ned(sLLA(1), sLLA(2)) * tECI2ECEF_mat;
        drNED = T * (pos - sPosECI);
        N = drNED(1); E = drNED(2);
        denom = max(N^2 + E^2, 1e-6);
        H(k, 1:3) = (-E * T(1,:) + N * T(2,:)) / denom;
    end
end