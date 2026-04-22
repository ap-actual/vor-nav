function plotStateEstimate(utc, x, P, trajImuData)
    % 1. Format time input
    utcSeconds = seconds(utc - utc(1));
    
    % 2. Define downsampling parameters
    FRAME_SKIP = 10;
    FRAME_END = numel(utcSeconds); 
    secIdx = 1:FRAME_SKIP:FRAME_END; 
    
    % 3. Downsample time and state data
    utcDownsampled = utc(secIdx);
    utcSecondsDownsampled = utcSeconds(secIdx);
    x_down = x(:, secIdx);
    P_down = P(:, :, secIdx);
    
    % 4. Downsample truth data (LLA)
    lat_down = trajImuData.lat(secIdx);
    lon_down = trajImuData.lon(secIdx);
    alt_down = trajImuData.alt(secIdx);
    
    % 5. Convert downsampled LLA to ECI
    % Pre-allocate to prevent the "zero-filling" issue
    numSamples = numel(secIdx);
    lla_truth = zeros(3, numSamples);
    
    for i = 1:numSamples
        % Convert each point using the corresponding timestamp
        lla_truth(:,i) = lla2eci([lat_down(i), lon_down(i), alt_down(i)], ...
                                 datevec(utcDownsampled(i)))';
    end
    
    % --- Figure 1: State vs Truth ---
    figure('Name', 'State Comparison (ECI)', 'NumberTitle', 'off');
    titles = ["ECI X", "ECI Y", "ECI Z"];
    for axisIdx = 1:3
        subplot(3,1,axisIdx); 
        plot(utcSecondsDownsampled, lla_truth(axisIdx,:), 'k', 'LineWidth', 1.5); hold on;
        plot(utcSecondsDownsampled, x_down(axisIdx,:), 'r--'); 
        legend("Truth", "Estimate"); 
        title(titles(axisIdx));
        ylabel('Meters');
        grid on;
    end
    xlabel('Time (seconds)');

    % --- Figure 2: Error and Covariance (Sigma Envelopes) ---
    figure('Name', 'Estimation Error (3-Sigma)', 'NumberTitle', 'off');
    for axisIdx = 1:3
        subplot(3,1,axisIdx);
        
        % Calculate error
        error_val = x_down(axisIdx,:) - lla_truth(axisIdx,:);
        
        % Calculate 1-sigma bounds from the diagonal of covariance matrix P
        sigma = sqrt(squeeze(P_down(axisIdx, axisIdx, :)))';
        
        % Plot error and 3-sigma bounds
        plot(utcSecondsDownsampled, error_val, 'b'); hold on;
        plot(utcSecondsDownsampled, sigma, 'r--'); 
        plot(utcSecondsDownsampled, -sigma, 'r--'); 
        
        title(['Error in ', titles(axisIdx)]);
        legend("Error", "1\sigma Bound");
        ylabel('Error (m)');
        grid on;
    end
    xlabel('Time (seconds)');
end