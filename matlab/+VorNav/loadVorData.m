function navAids = loadVorData(navAidFile)
arguments
    navAidFile {mustBeFile}
end

navAidTable = readtable(navAidFile);

% Filter for VOR classes
isLow = contains(navAidTable.CLASS_TXT, 'L-VOR');
isHigh = contains(navAidTable.CLASS_TXT, 'H-VOR');

% Convert to struct for faster indexing in loops
low = table2struct(navAidTable(isLow, {'X', 'Y', 'CLASS_TXT', 'IDENT'}));
high = table2struct(navAidTable(isHigh, {'X', 'Y', 'CLASS_TXT', 'IDENT'}));

% Standardize field names
low = arrayfun(@(s) struct('x', s.X, 'y', s.Y, 'class', s.CLASS_TXT, 'ident', s.IDENT), low);
high = arrayfun(@(s) struct('x', s.X, 'y', s.Y, 'class', s.CLASS_TXT, 'ident', s.IDENT), high);

navAids.low = low;
navAids.high = high;

fprintf('NAVAID Data Loaded: %d Low, %d High VORs.\n', length(navAids.low), length(navAids.high));

end