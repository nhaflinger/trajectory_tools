function downloadCoastlines()
%DOWNLOADCOASTLINES  One-time download of Natural Earth coastline data.
%
%   downloadCoastlines()
%
%   Fetches the Natural Earth 110m coastline GeoJSON from GitHub and saves
%   (coastlon, coastlat) as coastlines_cache.mat next to this file.
%   Handles both LineString and MultiLineString geometries, and both struct-
%   array and cell-array feature lists as returned by jsondecode.
%
%   Run once before using plotGroundTrack(), or it will be called automatically.

savePath = fullfile(fileparts(mfilename('fullpath')), 'coastlines_cache.mat');

url = ['https://raw.githubusercontent.com/nvkelso/' ...
       'natural-earth-vector/master/geojson/ne_110m_coastline.geojson'];

fprintf('Downloading Natural Earth 110m coastlines...\n');
try
    raw = webread(url, weboptions('Timeout', 30, 'ContentType', 'json'));
catch ME
    error('downloadCoastlines: download failed — %s\nCheck internet connection.', ME.message);
end

if ~isstruct(raw) || ~isfield(raw, 'features')
    error('downloadCoastlines: unexpected response format.');
end

% features may be a struct array or a cell array depending on field uniformity
features = raw.features;
if isstruct(features)
    nf = numel(features);
    getFeature = @(k) features(k);
else
    nf = numel(features);
    getFeature = @(k) features{k};
end

coastlon = [];
coastlat = [];

for k = 1:nf
    feat = getFeature(k);
    if ~isstruct(feat) || ~isfield(feat, 'geometry'), continue; end
    geom = feat.geometry;
    if isempty(geom) || ~isstruct(geom) || ~isfield(geom, 'type'), continue; end

    switch lower(geom.type)
        case 'linestring'
            [lo, la] = coordsToVec(geom.coordinates);
            coastlon = [coastlon; lo; NaN]; %#ok<AGROW>
            coastlat = [coastlat; la; NaN]; %#ok<AGROW>

        case 'multilinestring'
            rings = geom.coordinates;
            if ~iscell(rings), rings = {rings}; end
            for j = 1:numel(rings)
                [lo, la] = coordsToVec(rings{j});
                coastlon = [coastlon; lo; NaN]; %#ok<AGROW>
                coastlat = [coastlat; la; NaN]; %#ok<AGROW>
            end
    end
end

save(savePath, 'coastlon', 'coastlat');
fprintf('Saved %d coordinate points to:\n  %s\n', sum(~isnan(coastlon)), savePath);
end

%% ── Helper: convert jsondecode coordinate array to column vectors ───────────
function [lo, la] = coordsToVec(coords)
% coords is either:
%   - Nx2 double matrix        (jsondecode uniform array)
%   - Nx1 cell of 1x2 doubles  (jsondecode non-uniform array)
    if iscell(coords)
        coords = cell2mat(cellfun(@(c) c(:)', coords, 'UniformOutput', false));
    end
    lo = coords(:, 1);
    la = coords(:, 2);
end
