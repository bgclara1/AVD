function plotOmegaFrame(h, b, c, s, t)
    figure;
    hold on;
    
    % Define section coordinates
    % Webs (Vertical sections)
    webX = [-c/2, -c/2, -c/2+s, -c/2+s, c/2-s, c/2-s, c/2, c/2];
    webY = [0, h, h, 0, 0, h, h, 0];
    
    % Top flange (Horizontal section)
    topFlangeX = [-c/2, c/2, c/2, -c/2];
    topFlangeY = [h, h, h+t, h+t];

    % Bottom flanges (Horizontal extensions)
    bottomFlangeX = [-b/2, -c/2, -c/2, -b/2, b/2, c/2, c/2, b/2];
    bottomFlangeY = [0, 0, -t, -t, -t, -t, 0, 0];

    % Plot each section
    fill(webX, webY, 'b', 'FaceAlpha', 0.5); % Webs
    fill(topFlangeX, topFlangeY, 'b', 'FaceAlpha', 0.5); % Top flange
    fill(bottomFlangeX, bottomFlangeY, 'b', 'FaceAlpha', 0.5); % Bottom flanges

    % Formatting
    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Omega Section Geometry');
    axis equal;
    grid on;
    hold off;
end
