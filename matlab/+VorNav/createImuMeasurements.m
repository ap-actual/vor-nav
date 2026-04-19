function imuMeasured = createImuMeasurements(flightData, flightId)

    thisName = flightId;

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
    imuTruth.(thisName) = VorNav.buildImuFromFlightData( ...
        lat, lon, velocity, heading, vertrate, geoalt, time);

    % build imuSpec
    imuSpec.(thisName) = VorNav.buildImuSpec(); 
    
    % create imuMeasurements
    imuMeasured = VorNav.addImuErrors(imuTruth.(thisName),imuSpec.(thisName));

    % normalize IMU measurement times
    imuMeasured.time = imuMeasured.time - imuMeasured.time(1);

end
