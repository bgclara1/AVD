function plotIStringer(h, l, t)
    % Define coordinates for the I-stringer components
    
    % Web (vertical section)
    web_x = [-t/2, t/2, t/2, -t/2, -t/2];
    web_y = [-h/2, -h/2, h/2, h/2, -h/2];

    % Top flange
    top_x = [-l/2, l/2, l/2, -l/2, -l/2];
    top_y = [h/2, h/2, h/2 + t, h/2 + t, h/2];

    % Bottom flange
    bottom_x = [-l/2, l/2, l/2, -l/2, -l/2];
    bottom_y = [-h/2, -h/2, -h/2 - t, -h/2 - t, -h/2];

    figure;
    hold on;
    fill(web_x, web_y, 'r', 'FaceAlpha', 0.7, 'EdgeColor', 'k'); % Web in red
    fill(top_x, top_y, 'b', 'FaceAlpha', 0.7, 'EdgeColor', 'k'); % Top flange in blue
    fill(bottom_x, bottom_y, 'b', 'FaceAlpha', 0.7, 'EdgeColor', 'k'); % Bottom flange in blue
    
    axis equal;
    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Optimized I-Stringer Cross-Section');
    grid on;
end
