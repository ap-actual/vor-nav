function imu = addImuErrors(imu, imuSpec)

    %% accel
    accelTruth = imu.specificForceBody;
    N = size(accelTruth,1);

    imu.specificForceBodyTruth = accelTruth;

    % zero g level
    accelBias = imuSpec.accel.zeroGLevel * randn(1,3);

    % Linear acceleration sensitivity 2% 3 sigma
    accelScaleSigma = imuSpec.accel.linearAccelerationSensitivityPct / 3;
    accelScale = 1 + accelScaleSigma * randn(1,3);
    accelS = diag(accelScale);

    % Cross axis sensitivity
    accelM = eye(3) + imuSpec.accel.crossAxisSensitivity * [ ...
        0 1 1
        1 0 1
        1 1 0 ];

    % apply errors
    accelMeas = (accelM * (accelS * accelTruth.')).';
    accelMeas = accelMeas + accelBias;

    % saturate
    accelMeas = min(max(accelMeas, -imuSpec.accel.fullScale), imuSpec.accel.fullScale);

    % quantize using Linear acceleration sensitivity
    accelMeas = round(accelMeas / imuSpec.accel.linearAccelerationSensitivity) ...
              * imuSpec.accel.linearAccelerationSensitivity;

    imu.specificForceBodyMeas = accelMeas;

    %% gyro
    gyroTruth = imu.gyroRates;

    imu.gyroRatesTruth = gyroTruth;

    % Angular rate zero-rate level
    gyroBias = imuSpec.gyro.angularRateZeroRateLevel * randn(1,3);

    % Angular rate sensitivity, 2% 3 sigma
    gyroScaleSigma = imuSpec.gyro.angularRateSensitivityPct / 3;
    gyroScale = 1 + gyroScaleSigma * randn(1,3);
    gyroS = diag(gyroScale);

    % Cross axis sensitivity
    gyroM = eye(3) + imuSpec.gyro.crossAxisSensitivity * [ ...
        0 1 1
        1 0 1
        1 1 0 ];

    % Rate noise density 
    gyroNoiseSigma = imuSpec.gyro.rateNoiseDensity * sqrt(imuSpec.fs/2);
    gyroNoise = gyroNoiseSigma * randn(N,3);

    % apply errors
    gyroMeas = (gyroM * (gyroS * gyroTruth.')).';
    gyroMeas = gyroMeas + gyroBias + gyroNoise;

    % saturate
    gyroMeas = min(max(gyroMeas, -imuSpec.gyro.fullScale), imuSpec.gyro.fullScale);

    % quantize using Angular rate sensitivity
    gyroMeas = round(gyroMeas / imuSpec.gyro.angularRateSensitivity) ...
             * imuSpec.gyro.angularRateSensitivity;

    imu.gyroRatesMeas = gyroMeas;
end