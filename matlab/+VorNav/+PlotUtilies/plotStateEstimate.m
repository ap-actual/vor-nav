function plotStateEstimate(utc, x, P, trajImuData)

    % format time input
    utcSeconds = seconds(utc - utc(1));

    % downsample state data
    tol = 1e-9;
    secIdx = mod(utcSeconds, 2) < tol | mod(utcSeconds, 2) > (2 - tol);
    % secIdx = true(1,numel(utcSeconds));

    % utc time
    utcDownsampled = utc(secIdx);
    utcSecondsDownsampled = utcSeconds(secIdx);
    
    % downsample state data
    x = x(:,secIdx);
    P = P(:,:,secIdx);
    
    % downsample truth data
    lat = trajImuData.lat(secIdx);
    lon = trajImuData.lon(secIdx);
    alt = trajImuData.alt(secIdx);

    % convert downsampled lla to eci
    for i = 1:numel(lat)
        lla_truth(:,i) = lla2eci([lat(i),lon(i), alt(i)], datevec(utcDownsampled(i)))';
    end

    
    figure();
    subplot(3,1,1); 
    plot(utcSecondsDownsampled, lla_truth(1,:)); hold on;
    plot(utcSecondsDownsampled, x(1,:)); legend("truth", "estimate"); title("ECI X")

    subplot(3,1,2); 
    plot(utcSecondsDownsampled, lla_truth(2,:)); hold on;
    plot(utcSecondsDownsampled, x(2,:)); legend("truth", "estimate"); title("ECI Y")

    subplot(3,1,3); 
    plot(utcSecondsDownsampled, lla_truth(3,:)); hold on;
    plot(utcSecondsDownsampled, x(3,:)); legend("truth", "estimate"); title("ECI Z")

    figure();
    plot(utcSecondsDownsampled, x(1,:)-lla_truth(1,:)); hold on;
    plot(utcSecondsDownsampled, sqrt(squeeze(P(1,1,:))), 'r--'); 
    plot(utcSecondsDownsampled, -sqrt(squeeze(P(1,1,:))), 'r--'); 

end