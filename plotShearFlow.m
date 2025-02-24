function [theta, PlotSFlow, circle, skinThickness,totalSFlow] = plotShearFlow(shearYieldStressSkin)


    n = 3.75;
    W = 330800 * 9.81;
    Q = n * W;    % Radial load
    r = 3.195;     % Fuselage radius
    P = 0;        % Tangential load
    T = 0;        % Torque

    phi(1) = 180;
    for i = 2:37
        phi(i) = phi(i-1) - 10;
        if phi(i-1) == 0
            phi(i) = 350;
        end
    end
    
    
    q = (Q * sin(deg2rad(phi))) / (pi * r); 
    p = (P * cos(deg2rad(phi))) / (pi * r); 
    t = (T + P * r) / (2 * pi * r^2); 
    
    totalSFlow = q + p + t;
    
    theta = deg2rad(phi); 
    
    circle = ones(length(phi))*r*1e6;
    PlotSFlow = [];
    for i = 1:length(circle)
        if totalSFlow(i)+circle(i)<circle(i)
            PlotSFlow(i) = abs(totalSFlow(i)) + circle(i);
        else
            PlotSFlow(i) = (totalSFlow(i)) + circle(i);
        end 
    end
    
    % 
    % figure
    % polarplot(theta,PlotSFlow,'b-','LineWidth', 3)
    % hold on;
    % polarplot(theta,circle,'r-', 'LineWidth', 2)
    % ax = gca; 
    % ax.ThetaZeroLocation = 'top';
    % ax.ThetaDir = 'clockwise';
    % title('Shear Flow Distribution Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');
    % grid on;
    % legend('Shear Flow', 'Fuselage', 'Location', 'best');
    % 
    
    
    maxTotSF = max(abs(totalSFlow)); 
    skinThickness = maxTotSF/shearYieldStressSkin;

end
