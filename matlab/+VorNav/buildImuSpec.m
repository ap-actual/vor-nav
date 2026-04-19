function imuSpec = buildImuSpec()
% builds out the imuSpec for gyro and accel errors 
% Output:
%   imuSpec.accel
%   imuSpec.gyro

% ignore temp specs because no temperature data
    g0 = 9.81;
    fs = 10;
    imuSpec.fs = fs;
    imuSpec.g0 = g0;


    %% Accel
    imuSpec.accel.fullScale = 4 * g0;

    % Linear acceleration sensitivity 0.122 mg/LSB
    imuSpec.accel.linearAccelerationSensitivity = 0.122e-3 * g0; % m/s^2/LSB
    imuSpec.accel.linearAccelerationSensitivityPct = 0.02;  % ±2%

    % Zero g level 10 mg 
    imuSpec.accel.zeroGLevel = 10e-3 * g0; % m/s^2

    % Cross axis sensitivity ±0.5%
    imuSpec.accel.crossAxisSensitivity = 0.5 / 100;

    %% Gyro
    imuSpec.gyro.fullScale = deg2rad(500);

    % Angular rate sensitivity 17.50 mdps/LSB (500 dps mode)
    imuSpec.gyro.angularRateSensitivity = deg2rad(17.50e-3);  % rad/s/LSB
    imuSpec.gyro.angularRateSensitivityPct = 0.02;

    % Angular rate zero-rate  1 dps 
    imuSpec.gyro.angularRateZeroRateLevel = deg2rad(1); % rad/s

    % Cross axis sensitivity 1%
    imuSpec.gyro.crossAxisSensitivity = 1/ 100;

    % Rate noise density 5 mdps/sqrt(Hz)
    imuSpec.gyro.rateNoiseDensity = deg2rad(5e-3);  % rad/s/sqrt(Hz)

    % ARW 0.21 deg/sqrt(hr)
    imuSpec.gyro.angularRandomWalk = deg2rad(0.21)/ 60;  % rad/sqrt(s)
end