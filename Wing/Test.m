% Load data from provided Catchpole.txt
fid = fopen('Catchpole.txt', 'r');
if fid == -1
    error('Failed to open Catchpole.txt. Ensure the file exists in the current directory.');
end
data = fscanf(fid, '%c');
fclose(fid);
%%

% Parse data into a structured form
sections = regexp(data, '\((.*?)\)', 'match');
curves = regexp(data, '\(.*?\)[^\(]*', 'match');

thickness_data = {};
for i = 1:length(curves)
    ratio = str2double(regexprep(sections{i}, '[()]', ''));
    if isnan(ratio)
        ratio = inf;
    end
    values = textscan(regexprep(curves{i}, '\(.*?\)', ''), '%f,%f', 'Delimiter', ',');
    area_points = values{1};
    stress_points = values{2};
    thickness_data{end+1} = struct('thickness', ratio, 'area', area_points, 'stress', stress_points);
end

%%

% Initialize output vectors
stress_output = zeros(size(thickness_vec));
matched_thickness = zeros(size(thickness_vec));

for idx = 1:length(thickness_vec)
    t = thickness_vec(idx);
    a = area_vec(idx);

    % Find the closest thickness ratio
    differences = arrayfun(@(x) abs(x.thickness - t), [thickness_data{:}]);
    [~, min_idx] = min(differences);

    if t > 2
        min_idx = find([thickness_data{:}].thickness == inf, 1);
    end

    selected_curve = thickness_data{min_idx};

    % Perform linear interpolation for the given area ratio
    stress_output(idx) = interp1(selected_curve.area, selected_curve.stress, a, 'linear', 'extrap');
    matched_thickness(idx) = selected_curve.thickness;
    end