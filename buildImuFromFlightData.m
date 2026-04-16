function imu = buildImuFromFlightData(lat, lon, speed, heading, ...
    vertRate, geoAlt, time)
% Build approximate truth IMU data for a large airliner from flight track
% data. 
% 1. Take geodetic data (LLA) from a flight and interpolate it to chosen
%    imu rate.
% 2. Use interpolated data to calculate velocity in NED Nav frame
% 3. Calculate acceleration in NED Nav frame
% 4. Use simplifying arliner assumptions to get body euler rates P,Q,R for
%    gyro rates.
% 5. Convert acceleration into body frame for specf measurements
%
% Inputs:
%   lat:      (Nx1) latitude  [deg]
%   lon       (Nx1) longitude, [deg]
%   speed     (Nx1) total speed, [KNOTS]
%   heading   (Nx1) heading, [deg] (clockwise from north)
%   vertRate  (Nx1) vertical rate, [m/s] (up positive)
%   geoAlt    (Nx1) geometric altitude, [m]
%   time      (Nx1) sample times [sec]
%
% Output struct imu fields:
%    imu.time:              (Mx1) interpolated time to specified rate [sec]
%    imu.specificForceBody: (Mx3) dVs sensed in body frame [m/s^2]
%    imu.specificForceNav:  (Mx3) dvs in the nav frame [m/s^2]
%    imu.gyroRates:         (Mx3) pqr sensed in body frame [rad/s]
%    imu.eulerAngles:       (Mx3) body attitude relative to nav NED as euler
%                                 angles [rads] 
%    imu.velocityNED:        (Mx3) velocity in the nav frame [m/s]
%    imu.accelNED:           (Mx3) accelerations in nav frame (includes grav) [m/s^2]
%    imu.lat:              (Mx1) latitude [degs] 
%    imu.lon            (Mx1) longitude [degs]
%    imu.alt            (Mx1) altitude [m]

 
    % make sure these all columns
    lat = lat(:); % degrees 
    lon = lon(:); % degrees 
    speed = speed(:); % meter per sec
    heading = heading(:); % degrees 
    vertRate = vertRate(:); % meter per sec
    geoAlt = geoAlt(:); % meter
    time = time(:); % seconds
  
    % arbitrarily chose to sample imu data to 10 Hz
    sampleImu = 10; % 10 hz 
    timeImu = (time(1):1/sampleImu:time(end)).';
    M = numel(timeImu);

    % interpolate to the imu rate we want 
    latInterp = interp1(time, lat, timeImu, 'pchip');
    lonInterp = interp1(time, lon, timeImu, 'pchip');
    speedInterp = interp1(time, speed, timeImu, 'pchip');
    headingInterp = interp1(time, unwrap(deg2rad(heading)), timeImu, 'pchip'); 
    vertRateInterp = interp1(time, vertRate, timeImu, 'pchip');
    geoAltInterp = interp1(time, geoAlt, timeImu, 'pchip');
    
    speedInterp = speedInterp * 0.514444; % convert from knots to m/s

    dt = 1 / sampleImu;

    % get wgs84 gravity at lat and alt
    g  = gravitywgs84(latInterp,geoAltInterp);

    % given velocity is speed magnitude 
    vh = sqrt(max(speedInterp.^2 - vertRateInterp.^2, 0));

    % convert to NED frame so we are in the nav frame
    Vn = vh .* cos(headingInterp);
    Ve = vh .* sin(headingInterp);
    Vd = -vertRateInterp;
    velNav = [Vn Ve Vd];

    % calcualte acceleration from v and dt
    An = gradient(Vn, dt);
    Ae = gradient(Ve, dt);
    Ad = gradient(Vd, dt);
    accelNav = [An Ae Ad];


    % Make estimates of flight path angle and attitude using the 
    % following large airliner assumptions:
    % yaw (psi) = heading
    % pitch (theta) ~ flight path angle
    % roll (phi) from no slip turn
    psi = headingInterp;
    psiDot = gradient(unwrap(psi), dt);

    gamma = atan2(vertRateInterp, max(vh, 1e-6));   % up-positive flight path angle
    theta = gamma;                               
   
    phi = atan2(vh .* psiDot, g);

    % smooth attitude
    phi = smoothdata(phi,'movmean', 5);
    theta = smoothdata(theta,'movmean', 5);

    % euler angle rates
    phiDot = gradient(phi, dt);
    thetaDot = gradient(theta, dt);
    psiDot = psiDot; % already calculated above for phi

    % calculate PQR
    p = zeros(M,1);
    q = zeros(M,1);
    r = zeros(M,1);

    for k = 1:M
        p(k) = phiDot(k) - sin(theta(k)) * psiDot(k);
        q(k) = cos(phi(k)) * thetaDot(k) + sin(phi(k)) * cos(theta(k)) * psiDot(k);
        r(k) = -sin(phi(k)) * thetaDot(k) + cos(phi(k)) * cos(theta(k)) * psiDot(k);
    end

    % put everything from body frame to NED
    gyroBody = [p q r];
    gNav = [zeros(M,1), zeros(M,1), g];

    %-----------------------------
    % specific force in body frame:
    % f_b = Cbn * (a_ned - g_ned)
    % where Cbn maps NED -> body
    %-----------------------------
    specificForce = zeros(M,3);
    specificForceNav = zeros(M,3);
    euler = [phi theta psi];

    for k = 1:M
        tNav2Body = angle2dcm(psi(k),theta(k),phi(k),'ZYX'); % 321 rotation from nav 2 body
        specfNav = accelNav(k,:)' - gNav(k,:)'; % specific force in the nav frame 
        specfBody = tNav2Body * specfNav;
        specificForce(k,:) = specfBody.';
        specificForceNav(k,:) = specfNav.';
    end

    %-----------------------------
    % optional light smoothing
    %-----------------------------
    specificForce(:,1) = smoothdata(specificForce(:,1), 'movmean', 3);
    specificForce(:,2) = smoothdata(specificForce(:,2), 'movmean', 3);
    specificForce(:,3) = smoothdata(specificForce(:,3), 'movmean', 3);

    gyro(:,1) = smoothdata(gyroBody(:,1), 'movmean', 3);
    gyro(:,2) = smoothdata(gyroBody(:,2), 'movmean', 3);
    gyro(:,3) = smoothdata(gyroBody(:,3), 'movmean', 3);

    %-----------------------------
    % output
    %-----------------------------
    imu.time = timeImu;
    imu.specificForceBody = specificForce; % dVs sensed in body frame
    imu.gyroRates = gyro; % pqr sensed in body frame
    imu.specificForceNav = specificForceNav; % dvs in the nav frame
    imu.eulerAngles = euler; % body attitude relative to nav NED as euler angles
    imu.velocityNED = velNav; % velocity in the nav frame
    imu.accelNED = accelNav; % accelerations in nav frame 
    imu.lat = latInterp;
    imu.lon = lonInterp;
    imu.alt = geoAltInterp;
end