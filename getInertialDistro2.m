function [InertialLoads2] = getInertialDistro2()
    lTot = 80;

    fuselagePoints = [];
    fuselageWeight = 29744*9.81;
    for i = 1:lTot
        fuselagePoints(i) = fuselageWeight/lTot;
        xDiscr(i) = i-1;
    end
    
    payloadPoints = [];
    payloadWeight = 51879*9.81;
    % 
    % for i = 1:lTot
    %     if i<=6
    %         payloadPoints(i) = 0;
    %     elseif i>=61 
    %         payloadPoints(i) = 0;
    %     else 
    %         payloadPoints(i) = payloadWeight/55;
    %     end
    % end
    % 
    
    furnishingPoints = [];
    furnishingWeight = 17227*9.81;
    
    % for i = 1:lTot
    %     if i<=6
    %         furnishingPoints(i) = 0;
    %     elseif i>=61 
    %         furnishingPoints(i) = 0;
    %     else 
    %         furnishingPoints(i) = furnishingWeight/55;
    %     end
    % end
    % 
    fuelPoints = [];
    fuelWeight= 21000*9.81;
    
    % for i = 1:lTot
    %     if i<=31
    %         fuelPoints(i) = 0;
    %     elseif i>=39 
    %         fuelPoints(i) = 0;
    %     else 
    %         fuelPoints(i) = fuelWeight/8;
    %     end
    % end

    % components = {
    %     "Horizontal Tailplane", 2998, 75;
    %     "Vertical Tailplane", 1552.8, 72;
    %     "Nose Landing Gear", 486.52, 2;
    %     "Fuel Systems", 936, 46;
    %     "Flight Controls", 927, 47;
    %     "Installed APU", 1496, 71;
    %     "Instruments", 796.77, 3;
    %     "Hydraulic System", 271.9, 46;
    %     "Electronic System", 1037, 37;
    %     "Avionics", 697.64, 3;
    %     "Air Conditioning", 3593, 33; 
    %     "Anti Icing System", 757.9, 42;
    %     "Handling Gear", 113.3685, 33;
    % };
    % 
    % for i = 1:size(components)
    %     weights(i) = components{i, 2} * 9.81; % Mass (kg) * 9.81 to get weight (N)
    %     positions(i) = components{i, 3};      % x_cg position (m)
    % end
    % 
    % [uniquePositions, ~, idx] = unique(positions); % Unique positions and their indices
    % aggregatedWeights = accumarray(idx, weights); % Sum weights for each unique position


    % 
    % total = [];
    % 
    % for i = 1:lTot
    %     total(i) = fuselagePoints(i) + payloadPoints(i) + furnishingPoints(i) + fuelPoints(i);
    % end
    
%    total = fuselagePoints;


    % 
   %  total(42) = total(42) + 24237 * 9.81;    % Wing
   % % total(33) = total(33) + 29744 * 9.81;    % Fuselage
   %  total(75) = total(75) + 2998 * 9.81;     % Horizontal Tailplane
   %  total(72) = total(72) + 1552.8 * 9.81;   % Vertical Tailplane
   %  total(41) = total(41) + 4618 * 9.81;     % Engine Nacelle
   %  total(41) = total(41) + 23239 * 9.81;    % Engines
   %  total(41) = total(41) + 259 * 9.81;      % Engine Pneumatic Starter
   %  total(41) = total(41) + 1104 * 9.81;     % Engine Controls
   %  total(43) = total(43) + 5123 * 9.81;     % Main Landing Gear
   %  total(2)  = total(2)  + 486.52 * 9.81;   % Nose Landing Gear
   %  total(46) = total(46) + 936 * 9.81;      % Fuel Systems
   %  total(47) = total(47) + 927 * 9.81;      % Flight Controls
   %  total(71) = total(71) + 1496 * 9.81;     % Installed APU
   %  total(3)  = total(3)  + 796.77 * 9.81;   % Instruments
   %  total(46) = total(46) + 271.9 * 9.81;    % Hydraulic System
   %  total(37) = total(37) + 1037 * 9.81;     % Electronic System
   %  total(3)  = total(3)  + 697.64 * 9.81;   % Avionics
   %  total(33) = total(33) + 17227 * 9.81;    % Furnishing
   %  total(33) = total(33) + 3593 * 9.81;     % Air Conditioning
   %  total(42) = total(42) + 757.9 * 9.81;    % Anti Icing System
   %  total(33) = total(33) + 113.3685 * 9.81; % Handling Gear
   %  total(33) = total(33) + 51879 * 9.81;    % Payload
   %  total(35) = total(35) + fuelWeight;

   % Component positions (sorted)
x = [2, 3, 33, 35, 37, 41, 42, 43, 46, 47, 71, 72, 75];

% Corresponding total component weights at positions (sorted)
y = [
    486.52 * 9.81;                             % Nose Landing Gear
    (796.77 + 697.64) * 9.81;                  % Instruments + Avionics
    (17227 + 3593 + 113.3685 + 51879) * 9.81;  % Furnishing + Air Cond. + Handling Gear + Payload
    156760 * 9.81;                             % Fuel (total)
    1037 * 9.81;                               % Electronic System
    (4618 + 23239 + 259 + 1104) * 9.81;        % Engine Nacelle + Engines + Pneumatic Starter + Controls
    (24237 + 757.9) * 9.81;                    % Wing + Anti Icing System
    5123 * 9.81;                               % Main Landing Gear
    (936 + 271.9) * 9.81;                      % Fuel Systems + Hydraulic System
    927 * 9.81;                                % Flight Controls
    1496 * 9.81;                               % Installed APU
    1552.8 * 9.81;                             % Vertical Tailplane
    2998 * 9.81                                % Horizontal Tailplane
];

% Plotting the discrete weights as vertical lines for clarity
    figure;
    bar(x, y*-1,'LineWidth', 1.5);
    xlabel('Fuselage Position (m)');
    ylabel('Component Weight (N)');
    title('Discrete Component Weight Distribution');
    grid on;


    % 
    % total = total*-1;
    % InertialLoads2 = total;

    % figure
    % bar(xDiscr, total, 'b-', 'LineWidth', 1.5);
    % xlabel('Fuselage Position (m)','FontWeight','bold');
    % ylabel('Load (N/m)','FontWeight','bold');
    % title('Inertial Loads');
    % grid minor
    % 

end