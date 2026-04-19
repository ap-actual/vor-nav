function velIntNED = integrateSpecificForceNav(imu)
    % g = 9.80665; using this gives big delta
    N = length(imu.time);
    g = gravitywgs84(imu.lat,imu.alt);

    % recover true acceleration
    aNav = imu.specificForceNav + [zeros(N,1) zeros(N,1) g];

    % integrate
    velIntNED = imu.velocityNED(1,:) + ...
                cumsum(aNav * mean(diff(imu.time)), 1);

end