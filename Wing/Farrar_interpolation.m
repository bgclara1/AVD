function [farrar_output, matched_efficiency] = Farrar_interpolation(thickness_vec, area_vec)
    fid = fopen('Farrar.txt', 'r');
    if fid == -1
        error('Failed to open Farrar.txt. Ensure the file exists in the current directory.');
    end
    data = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    data = data{1};

    efficiency_data = {};
    current_efficiency = [];
    thickness_points = [];
    area_points = [];

    for i = 1:length(data)
        line = strtrim(data{i});
        if startsWith(line, '(')
            if ~isempty(current_efficiency)
                efficiency_data{end+1} = struct('efficiency', current_efficiency, 'thickness', thickness_points, 'area', area_points);
            end
            current_efficiency = str2double(regexprep(line, '[()]', ''));
            thickness_points = [];
            area_points = [];
        elseif ~isempty(line)
            values = strsplit(line, ',');
            if length(values) == 2
                thickness_val = str2double(values{1});
                area_val = str2double(values{2});
                if isfinite(thickness_val) && isfinite(area_val)
                    thickness_points(end+1) = thickness_val;
                    area_points(end+1) = area_val;
                end
            end
        end
    end
    efficiency_data{end+1} = struct('efficiency', current_efficiency, 'thickness', thickness_points, 'area', area_points);

    farrar_output = zeros(size(thickness_vec));
    matched_efficiency = zeros(size(thickness_vec));

    for idx = 1:length(thickness_vec)
        t = thickness_vec(idx);
        a = area_vec(idx);
        distances = arrayfun(@(x) min(sqrt((x.thickness - t).^2 + (x.area - a).^2)), [efficiency_data{:}]);
        [~, min_idx] = min(distances);
        selected_curve = efficiency_data{min_idx};
        farrar_output(idx) = selected_curve.efficiency;
        matched_efficiency(idx) = selected_curve.efficiency;
    end
end

