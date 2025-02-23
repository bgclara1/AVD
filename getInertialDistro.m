
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           INERTIAL LOAD DISCRETISATION AND PLOTTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [xDiscr,total] = getInertialDistro()
    lTot = 80;
    DiscInterval = 5;
    fuselageDisc = 0:DiscInterval:lTot;
    cabinLow = 6.23;
    cabinHigh = 60.48;
    numberOfPoints = lTot/DiscInterval;
    
    fuselagePoints = [];
    fuselageWeight = 29744*9.81;
    for i = 1:lTot
        fuselagePoints(i) = fuselageWeight/lTot;
        xDiscr(i) = i-1;
    end
    
    payloadPoints = [];
    payloadWeight = 51879*9.81;
    
    for i = 1:lTot
        if i<=6
            payloadPoints(i) = 0;
        elseif i>=61 
            payloadPoints(i) = 0;
        else 
            payloadPoints(i) = payloadWeight/lTot;
        end
    end
    
    
    furnishingPoints = [];
    furnishingWeight = 17227*9.81;
    
    for i = 1:lTot
        if i<=6
            furnishingPoints(i) = 0;
        elseif i>=61 
            furnishingPoints(i) = 0;
        else 
            furnishingPoints(i) = furnishingWeight/lTot;
        end
    end
    
    fuelPoints = [];
    fuelWeight= 21000*9.81;
    
    for i = 1:lTot
        if i<=31
            fuelPoints(i) = 0;
        elseif i>=39 
            fuelPoints(i) = 0;
        else 
            fuelPoints(i) = fuelWeight/lTot;
        end
    end
    
    components = {
        "Horizontal Tailplane", 2998, 75;
        "Vertical Tailplane", 1552.8, 72;
        "Nose Landing Gear", 486.52, 2;
        "Fuel Systems", 936, 46;
        "Flight Controls", 927, 47;
        "Installed APU", 1496, 71;
        "Instruments", 796.77, 3;
        "Hydraulic System", 271.9, 46;
        "Electronic System", 1037, 37;
        "Avionics", 697.64, 3;
        "Air Conditioning", 3593, 33; 
        "Anti Icing System", 757.9, 42;
        "Handling Gear", 113.3685, 33;
    };
    
    for i = 1:size(components)
        weights(i) = components{i, 2} * 9.81; % Mass (kg) * 9.81 to get weight (N)
        positions(i) = components{i, 3};      % x_cg position (m)
    end
    
    [uniquePositions, ~, idx] = unique(positions); % Unique positions and their indices
    aggregatedWeights = accumarray(idx, weights); % Sum weights for each unique position
    
    
    
    total = [];

    for i = 1:lTot
        total(i) = fuselagePoints(i) + payloadPoints(i) + furnishingPoints(i) + fuelPoints(i);
    end
    % 
    % for i = 1:lTot
    %     total(i) = fuselagePoints(i);
    % end
    % 
    % 
    % 
    total(75) = total(75) + 2998*9.81; %HT
    total(72) = total(72) + 1552.8*9.81; %VT
    total(2) = total(2) + 486.52*9.81; %nose landing gear
    total(46) = total(46) + 936*9.81; %fuel sus
    total(47) = total(47) + 927*9.81; %flight controls
    total(3) = total(3) + 796.77*9.81; %instruments
    total(71) = total(71) + 1496*9.81; %APU
    total(37) = total(37) + 1037*9.81; %Elec 
    total(42) = total(42)+ 757.9*9.81; %anti icing
    
    total = total*-1;
    total = total*9.81;

    figure
    plot(xDiscr, total, 'b-', 'LineWidth', 1.5);
    
    xlabel('Fuselage Position (m)','FontWeight','bold');
    ylabel('Load (N/m)','FontWeight','bold');
    title('Inertial Loads');
    grid minor
 

end

