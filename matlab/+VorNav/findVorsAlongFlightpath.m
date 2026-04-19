function findVorsAlongFlightpath(trajImuData, navAids)

% get all 1-second time steps
tol = 1e-9;
secIdx = mod(trajImuData.time, 2) < tol | mod(trajImuData.time, 2) > (2 - tol);

t   = trajImuData.time(secIdx);
lat = trajImuData.lat(secIdx);
lon = trajImuData.lon(secIdx);
alt = trajImuData.alt(secIdx);


end