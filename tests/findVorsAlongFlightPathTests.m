
% load vor data
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% load flight data

[trajImuData, imuTruth] = createImuMeasurements(flightData, 'icao_abe79d_F1');

% 
VorNav.findVorsAlongFlightpath(trajImuData, navAids);