
% load vor data
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% load flight data
load('../dataSets/filtered_flights.mat');
trajImuData = VorNav.createImuMeasurements(flightData, 'icao_abe79d_F1');

% find vors along flight path
visibleVorIdents = VorNav.findVorsAlongFlightpath(trajImuData, navAids, 'plotFlag', true);

