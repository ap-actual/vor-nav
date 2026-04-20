function visibleVorIdents = findVorsAlongFlightpath(trajImuData, navAids, args)
arguments
    trajImuData struct
    navAids
    args.plotFlag = false
end

% get all 1-second time steps
tol = 1e-9;
secIdx = mod(trajImuData.time, 2) < tol | mod(trajImuData.time, 2) > (2 - tol);

% extract low-rate data
t   = trajImuData.time(secIdx);
lat = trajImuData.lat(secIdx);
lon = trajImuData.lon(secIdx);
alt = trajImuData.alt(secIdx);

visibleVorIdents = string.empty;

for i = 1:numel(t)
    
    vorData = VorNav.vorMeas([lat(i), lon(i), alt(i)], "lla", navAids);
    identsTemp = {vorData.ident};

    uniqueNew = setdiff(string(identsTemp), visibleVorIdents, 'stable');
    visibleVorIdents = [visibleVorIdents, uniqueNew];

end

visibleVorIdents = visibleVorIdents';


% plot if requested
if args.plotFlag
    figure();
    geoscatter(lat, lon); hold on;
    
    for i = 1:numel(navAids.high)
        % plot high VORs
        if ismember(string(navAids.high(i).ident), visibleVorIdents)
            geoscatter(navAids.high(i).y, navAids.high(i).x, 'or');
        end
    end
    
    for i = 1:numel(navAids.low)
        % plot low VORs
        if ismember(string(navAids.low(i).ident), visibleVorIdents)
            geoscatter(navAids.low(i).y, navAids.low(i).x, 'ob');
        end
    end
    geobasemap topographic
end



end