function [stress_output, matched_thickness] = catchpole_interpolation(thickness_vec, area_vec)
    fid = fopen('Catchpole.txt', 'r');
    if fid == -1
        error('Failed to open Catchpole.txt. Ensure the file exists in the current directory.');
    end
    data = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    data = data{1};

    thickness_data = {};
    current_thickness = [];
    area_points = [];
    stress_points = [];

    for i = 1:length(data)
        line = strtrim(data{i});
        if startsWith(line, '(')
            if ~isempty(current_thickness)
                thickness_data{end+1} = struct('thickness', current_thickness, 'area', area_points, 'stress', stress_points);
            end
            current_thickness = str2double(regexprep(line, '[()]', ''));
            if isnan(current_thickness)
                current_thickness = inf;
            end
            area_points = [];
            stress_points = [];
        elseif ~isempty(line)
            values = strsplit(line, ',');
            if length(values) == 2
                area_val = str2double(values{1});
                stress_val = str2double(values{2});
                if isfinite(area_val) && isfinite(stress_val)
                    area_points(end+1) = area_val;
                    stress_points(end+1) = stress_val;
                end
            end
        end
    end
    thickness_data{end+1} = struct('thickness', current_thickness, 'area', area_points, 'stress', stress_points);

    stress_output = zeros(size(thickness_vec));
    matched_thickness = zeros(size(thickness_vec));

    for idx = 1:length(thickness_vec)
        t = thickness_vec(idx);
        a = area_vec(idx);
        differences = arrayfun(@(x) abs(x.thickness - t), [thickness_data{:}]);
        [~, min_idx] = min(differences);
        if t > 2
            for k = 1:length(thickness_data)
                if isinf(thickness_data{k}.thickness)
                    min_idx = k;
                    break;
                end
            end
        end
        selected_curve = thickness_data{min_idx};
        if all(isfinite(selected_curve.area)) && all(isfinite(selected_curve.stress))
            stress_output(idx) = interp1(selected_curve.area, selected_curve.stress, a, 'linear', 'extrap');
        else
            stress_output(idx) = NaN;
        end
        matched_thickness(idx) = selected_curve.thickness;
    end
end
