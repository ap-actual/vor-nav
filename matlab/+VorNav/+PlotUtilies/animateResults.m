function animateResults(matFile)

%% 1. Load Saved Data
load(matFile, 'results', 'trajImuData', 'navAids');
fprintf('Data loaded from kf_results.mat\n');

%% 2. GIF Export Settings
GIF_FILENAME   = 'kf_output.gif';
GIF_FRAME_SKIP = 10;
GIF_DELAY      = 0.1;
GIF_RESOLUTION = '-r72';
GIF_COLORS     = 128;
JUMP_FRAMES    = 20;

%% 3. Pre-compute full truth LLA trajectory
nTimes   = numel(results.time);
truthLLA = [trajImuData.lat(:), trajImuData.lon(:)];  % Nx2, used for trail

%% 4. Animation Loop
fprintf('--- Starting Geoplot Animation ---\n');
gif_frame_count = 0;

for i = 1:JUMP_FRAMES:nTimes-1

    currentUTC = results.utc(i);

    % --- Determine visible VORs at this timestep ---
    truthLla_i = [trajImuData.lat(i), trajImuData.lon(i), trajImuData.alt(i)];
    vorMeasData = VorNav.vorMeas(truthLla_i, "lla", navAids);
    nVor        = numel(vorMeasData);

    % --- Convert KF estimate ECI -> LLA ---
    estLLA = eci2lla(results.x(1:3, i+1)', datevec(currentUTC));

    % --- Reconstruct nUsed from saved results ---
    nUsed = results.nUsed(i);

    % ------ GEOPLOT UPDATE -------------------------------------------
    figure(1); clf;

    % -- Plot full truth trail (faint) --
    geoplot(truthLLA(1:i,1), truthLLA(1:i,2), 'r:', ...
        'LineWidth', 0.8, 'DisplayName', 'Truth trail');
    hold on;

    % -- Active VOR stations --
    if nVor > 0
        activeLats = [vorMeasData.lat];
        activeLons = [vorMeasData.lon];
        geoplot(activeLats, activeLons, 'b^', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'b', 'DisplayName', 'Active VORs');

        for j = 1:nVor
            text(activeLats(j), activeLons(j), ...
                sprintf('  %s', vorMeasData(j).ident), ...
                'FontSize', 8, 'Color', 'blue');
        end
    end

    % -- KF Estimate --
    geoplot(estLLA(1), estLLA(2), 'g*', 'MarkerSize', 14, ...
        'MarkerFaceColor', 'g', 'DisplayName', 'KF Estimate');

    % -- Truth current position --
    geoplot(truthLla_i(1), truthLla_i(2), 'r+', 'MarkerSize', 12, ...
        'LineWidth', 2, 'DisplayName', 'Truth');

    % -- VOR Bearing Cones --
    if nVor > 0
        CONE_HALF_ANGLE_DEG = 2;

        for j = 1:nVor
            vorLat = vorMeasData(j).lat;
            vorLon = vorMeasData(j).lon;

            dLat = estLLA(1) - vorLat;
            dLon = estLLA(2) - vorLon;
            bearingToEst = atan2(dLon, dLat);

            CONE_LENGTH_DEG = sqrt(dLat^2 + dLon^2);

            leftBearing  = bearingToEst - deg2rad(CONE_HALF_ANGLE_DEG);
            rightBearing = bearingToEst + deg2rad(CONE_HALF_ANGLE_DEG);

            arcAngles = linspace(leftBearing, rightBearing, 30);
            arcLats   = vorLat + CONE_LENGTH_DEG * cos(arcAngles);
            arcLons   = vorLon + CONE_LENGTH_DEG * sin(arcAngles);

            coneLats = [vorLat, arcLats, vorLat];
            coneLons = [vorLon, arcLons, vorLon];

            geoplot(coneLats, coneLons, 'c-', 'LineWidth', 1, 'DisplayName', '');
        end

        geoplot(NaN, NaN, 'c-', 'LineWidth', 1.5, 'DisplayName', '2° VOR cone');
    end

    % -- 3-sigma Covariance Ellipse --
    P_pos_ECI = results.P(1:3, 1:3, i+1);

    tECI2ECEF = dcmeci2ecef('IAU-2000/2006', datevec(currentUTC));
    tECEF2NED = dcmecef2ned(estLLA(1), estLLA(2));
    tECI2NED  = tECEF2NED * tECI2ECEF;

    P_NED = tECI2NED * P_pos_ECI * tECI2NED';

    METERS_PER_DEG = 111139;
    P_latlon = P_NED(1:2, 1:2) / METERS_PER_DEG^2;

    [V, D]      = eig(P_latlon);
    angles      = linspace(0, 2*pi, 100);
    nSigma      = 3;
    unit_circle = [cos(angles); sin(angles)];
    ellipse_pts  = V * (nSigma * sqrt(abs(D))) * unit_circle;
    ellipse_lats = estLLA(1) + ellipse_pts(1, :);
    ellipse_lons = estLLA(2) + ellipse_pts(2, :);

    geoplot(ellipse_lats, ellipse_lons, 'y-', 'LineWidth', 1.5, ...
        'DisplayName', '3\sigma bounds');

    % -- Formatting --
    geobasemap('streets');
    legend('Location', 'best');
    title(sprintf('KF State  |  t = %.1f s  |  %d/%d bearings used', ...
        results.time(i), nUsed, nVor));
    drawnow;

    % -- GIF Frame Capture --
    gif_frame_count = gif_frame_count + 1;
    if mod(gif_frame_count, GIF_FRAME_SKIP) == 0
        frame = getframe(figure(1));
        im    = frame2im(frame);
        [imind, cm] = rgb2ind(im, GIF_COLORS);

        if gif_frame_count == GIF_FRAME_SKIP
            imwrite(imind, cm, GIF_FILENAME, 'gif', ...
                'Loopcount', Inf, ...
                'DelayTime', GIF_DELAY);
        else
            imwrite(imind, cm, GIF_FILENAME, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', GIF_DELAY);
        end
    end

end

fprintf('Animation complete. GIF saved to: %s\n', GIF_FILENAME);