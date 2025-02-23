function [x,y,numStringers] = stringerPlot(StringerSpacing)

    n = 3.75;
    W = 330800 * 9.81;
    r = 3.195;     % Fuselage radius

    % Given Data
    C = 2 * pi * r;
   % StringerSpacing = 0.5; % Change
    numStringers = floor(C / StringerSpacing);

    angle = linspace(0, 360, numStringers + 1);
    
    x = r * cosd(angle); 
    y = r * sind(angle); 
    
    % figure;
    % scatter(x, y, 'filled');
    % grid on;
    % axis equal; 
    % margin = 0.3 * max(abs([x, y])); 
    % xlim([-max(abs(x)) - margin, max(abs(x)) + margin]);
    % ylim([-max(abs(y)) - margin, max(abs(y)) + margin]);
    % xlabel('X-axis', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('Y-axis', 'FontSize', 12, 'FontWeight', 'bold');
    % title('Stringer Positions Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');

end