function [h, b, t, s] = lightFrameOptimisation(M_ult)
     
    % Material and structural properties
    E = 73*1e9;  
    D = 6.485;    
    %L = 0.516;   
    C_f = 6.25*1e-5; 
    % 
    % EI = C_f * M_ult * D^2 / L;
    % I_required = EI / E;

         % h,    b,     t,   s (spacing)
    lb = [0.02, 0.01, 0.001, 0.45];
    ub = [0.12, 0.5, 0.01, 0.6]; 

    hRange = linspace(lb(1), ub(1), 40);
    bRange = linspace(lb(2), ub(2), 40);
    tRange = linspace(lb(3), ub(3), 40);
    sRange = linspace(lb(4), ub(4), 40);

    success = [];

    fprintf('Total iterations: %d\n', length(hRange));

    for i = 1:length(hRange)
        disp(i)
        currentH = hRange(i);
        for j = 1:length(bRange)
            currentB = bRange(j);
            for k = 1:length(tRange)
                currentT = tRange(k);
                for l = 1:length(sRange)
                    currentS = sRange(l);

                   
                    if (currentH / currentB < 1.3) || (currentH / currentB > 5)
                        continue; 
                    end
                    if (currentH / currentT > 30)
                        continue;
                    end
                    if (currentB / currentT > 30)
                        continue;
                    end

                    EI = C_f * M_ult * D^2 / currentS;
                    I_required = EI / E;

                    currentIxx = (currentB * currentH^3) / 12 - ((currentB - currentT) * (currentH - 2*currentT)^3) / 12;

                    if currentIxx < I_required
                        continue;
                    end
                    
                    currentFrameArea = 2 * (currentH * currentT) + (currentB * currentT);
                 %   currentFrameArea = currentB*currentH-(currentH-2*currentT)*(currentB-currentT);
                    success = [success; i, j, k, l, currentFrameArea, currentIxx];

                end
            end
        end
    end

    % Extract successful values
    successfulH = hRange(success(:,1));
    successfulB = bRange(success(:,2));
    successfulT = tRange(success(:,3));
    successfulS = sRange(success(:,4));
    successfulArea = success(:,5);
    successfulIxx = success(:,6);

    % Optimisation objective: Minimise area while meeting stiffness
    objectiveValue = successfulArea'.*successfulS;
    [~, idxMin] = min(objectiveValue);

    % Plot results
    figure;
    scatter3(successfulH, successfulB, successfulT, 50, objectiveValue, 'o', ...
         'MarkerFaceColor', 'flat', 'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);
    xlabel('Frame Height (h)');
    ylabel('Frame Width (b)');
    zlabel('Frame Thickness (t)');
    colorbar;
    title('Light Frame C-Section Optimisation');

    hold on;
    plot3(successfulH(idxMin), successfulB(idxMin), successfulT(idxMin),...
        'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    legend('Design Points', 'Minimum Mass');

    % Return optimised dimensions
    h = successfulH(idxMin);
    b = successfulB(idxMin);
    t = successfulT(idxMin);
    s = successfulS(idxMin);

end
