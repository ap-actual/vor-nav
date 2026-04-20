function validVORs = vorMeas(pos, posType, navaids)
    % VORMEAS Returns VORs in range and bearing
    arguments
        pos (1,3) double
        posType {mustBeMember(posType, ["lla", "ecef"])}
        navaids struct
    end

    % Constants
    nm2m = 1852;
    validVORs = [];
    m2ft = 3.28084;

    % Normalize to row vector regardless of input shape
    pos = pos(:)';

    % 1. Coordinate Handling
    if posType == "ecef"
        myLatLon = ecef2lla(pos);
        myLat = myLatLon(1);
        myLon = myLatLon(2);
        myAlt = myLatLon(3) * m2ft;
    else
        myLat = pos(1);
        myLon = pos(2);
        myAlt = pos(3) * m2ft;
    end

    % 2. Calculation Handles
    haversine = @(lat1, lon1, lat2, lon2) ...
        2 * 6371000 * asin(sqrt(sin(deg2rad((lat2 - lat1)/2)).^2 + ...
        cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(deg2rad((lon2 - lon1)/2)).^2));

    % Calculates bearing from Point A to Point B
    calculateBearing = @(latA, lonA, latB, lonB) ...
        mod(rad2deg(atan2(sin(deg2rad(lonB-lonA)) * cos(deg2rad(latB)), ...
        cos(deg2rad(latA)) * sin(deg2rad(latB)) - ...
        sin(deg2rad(latA)) * cos(deg2rad(latB)) * cos(deg2rad(lonB-lonA)))), 360);

    % 3. Service Volume Logic
    highRanges = [1000, 14500, 40; 14500, 18000, 100; 18000, 45000, 130; 45000, 60000, 100];
    rangeNM = 0;
    for i = 1:size(highRanges,1)
        if myAlt >= highRanges(i,1) && myAlt <= highRanges(i,2)
            rangeNM = highRanges(i,3);
            break;
        end
    end

    % 4. Search VORs
    allTypes = [navaids.low; navaids.high];
    for i = 1:length(allTypes)

        % Explicitly map x/y to lat/lon for clarity
        vorLat = allTypes(i).y;
        vorLon = allTypes(i).x;

        isLow = contains(allTypes(i).class, 'L-VOR');
        currentMaxRange = rangeNM;
        if isLow, currentMaxRange = 40; end
        if isLow && myAlt >= 14500, continue; end

        dist = haversine(myLat, myLon, vorLat, vorLon);

        if dist <= currentMaxRange * nm2m
            res = allTypes(i);

            % Add explicit coordinate fields to the output struct
            res.lat = vorLat;
            res.lon = vorLon;
            res.distance_m = dist;

            % Bearing FROM VOR TO aircraft (true VOR radial)
            res.bearing_deg = calculateBearing(vorLat, vorLon, myLat, myLon);

            validVORs = [validVORs; res];
        end
    end
end