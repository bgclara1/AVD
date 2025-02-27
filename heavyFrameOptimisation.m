function [h, b, c,  t] = heavyFrameOptimisation(maxN, maxS, maxM)
     
    tensileYieldStressFrame = 310e6;
    shearYieldStressFrame = 250e6;

     % [Web height, Bottom flange width, Top flange width, Base thickness]
    lb = [0.001, 0.001, 0.001, 0.003]; 
    ub = [0.5, 0.5,0.5, 0.02];    

    hRange = linspace(lb(1), ub(1), 75);
    bRange = linspace(lb(2), ub(2), 75);
    cRange = linspace(lb(3), ub(3), 75);
    tRange = linspace(lb(4), ub(4), 75);

    success = [];
    for i = 1:length(hRange)
        currentH = hRange(i);
        for j = 1:length(bRange)
            currentB = bRange(j);
            for k = 1:length(cRange)
                currentC = cRange(k);
                for l = 1:length(tRange)
                    currentT = tRange(l);
                 % New constraints for currentC
                    if (currentC / currentH < 1.5) || (currentC / currentH > 4)
                        continue; % Ensure currentC is 1.5 to 4 times bigger than currentH
                    end
                    if (currentC / currentB < 1.5) || (currentC / currentB > 4)
                        continue; % Ensure currentC is 1.5 to 4 times bigger than currentB
                    end
                    if (currentC / currentT < 4) || (currentC / currentT > 15)
                        continue; % Ensure currentC is 1.5 to 4 times bigger than currentB
                    end
                                        % Compute cross-section properties
                    [currentFrameArea, currentFrameIxx] = OmegaFrameProps(currentH, currentB, currentC, currentT);

                    % Structural Compliance Checks
                    directStressCompliant = maxN <= tensileYieldStressFrame * currentFrameArea;
                    shearStressCompliant = maxS <= shearYieldStressFrame * currentFrameArea;
                    bendingStressCompliant = (maxM/currentFrameIxx)*(currentH/2) <= tensileYieldStressFrame ;

                    structurallyCompliant = all([directStressCompliant, shearStressCompliant, bendingStressCompliant]);

                    if structurallyCompliant
                        success = [success; i, j, k, l, currentFrameArea];
                    end
             
                end
            end
        end
    end

    % Extract successful values
    successfulH = hRange(success(:,1));
    successfulB = bRange(success(:,2));
    successfulC = cRange(success(:,3));
    successfulT = tRange(success(:,4));
    successfulArea = success(:,5);

    % **Objective function: Minimize cross-section area**
    objectiveValue = successfulArea;
    [~, idxMin] = min(objectiveValue);

    % **3D Scatter Plot of Optimized Omega Sections**
    figure;
    scatter3(successfulH, successfulB, successfulC,50, objectiveValue, 'filled');
    xlabel('Web Height (H)');
    ylabel('Bottom Flange Width (B)');
    zlabel('Top Flange Width (C)');
    colorbar;
    title('Heavy Frame Omega-Section Optimisation');

    hold on;
    plot3(successfulH(idxMin), successfulB(idxMin), successfulC(idxMin),...
        'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    legend('Design Points', 'Minimum Objective');

    % **Assign optimal values**
    h = successfulH(idxMin);
    b = successfulB(idxMin);
    c = successfulC(idxMin);
    t = successfulT(idxMin);

    % **Plot Omega Section**
   % plotOmegaFrame(h, b, c, s, t);
end
