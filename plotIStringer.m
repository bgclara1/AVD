function plotIStringer(h, l, t)

    web_x = [-t/2, t/2, t/2, -t/2, -t/2];
    web_y = [-h/2, -h/2, h/2, h/2, -h/2];


    top_x = [-l/2, l/2, l/2, -l/2, -l/2];
    top_y = [h/2, h/2, h/2 + t, h/2 + t, h/2];


    bottom_x = [-l/2, l/2, l/2, -l/2, -l/2];
    bottom_y = [-h/2, -h/2, -h/2 - t, -h/2 - t, -h/2];

    figure;
    hold on;
    fill(web_x, web_y, 'r', 'FaceAlpha', 0.7, 'EdgeColor', 'k');d
    fill(top_x, top_y, 'b', 'FaceAlpha', 0.7, 'EdgeColor', 'k'); 
    fill(bottom_x, bottom_y, 'b', 'FaceAlpha', 0.7, 'EdgeColor', 'k'); 
    
    axis equal;
    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Optimized I-Stringer Cross-Section');
    grid on;
end
