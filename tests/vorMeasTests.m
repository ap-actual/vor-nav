% this script tests that the vorMeas function works
navAids = VorNav.loadVorData('../dataSets/NAVAID_System.csv');

% test LLA
curLLA     = [40.5, -80.2, 10000];
resultsLla = VorNav.vorMeas(curLLA, "lla", navAids);

% test ECEF
curEcef     = lla2ecef(curLLA);
resultsEcef = VorNav.vorMeas(curEcef, "ecef", navAids);