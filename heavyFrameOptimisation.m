function [h, b, c,  t] = heavyFrameOptimisation(maxN, maxS, maxM)
     
    tensileYieldStressFrame = 413e6;
    shearYieldStressFrame = 250e6;

     % [Web height, Bottom flange width (2 short), Top flange width (long) , Base thickness]
    lb = [0.001, 0.001, 0.001, 0.001]; 
    ub = [0.3 0.3, 0.3, 0.1];    

    hRange = linspace(lb(1), ub(1), 30);
    bRange = linspace(lb(2), ub(2), 30);
    cRange = linspace(lb(3), ub(3), 30);
    tRange = linspace(lb(4), ub(4), 30);

    success = [];
    for i = 1:length(hRange)
        currentH = hRange(i);
        for j = 1:length(bRange)
            currentB = bRange(j);
            for k = 1:length(cRange)
                currentC = cRange(k);
                for l = 1:length(tRange)
                    currentT = tRange(l);
          
                    if (currentC / currentH < 1.5) || (currentC / currentH > 15)
                        continue; 
                    end
                    if (currentC / currentB < 1.2) || (currentC / currentB > 15)
                        continue; 
                    end
                    if (currentC / currentT < 3.5) || (currentC / currentT > 30)
                        continue; 
                    end
                               
                    [currentFrameArea, currentFrameIxx] = OmegaFrameProps(currentH, currentB, currentC, currentT);

                    directStressCompliant = maxN/currentFrameArea <= tensileYieldStressFrame ;
                    shearStressCompliant = maxS/currentFrameArea <= shearYieldStressFrame;
                    bendingStressCompliant = (maxM/currentFrameIxx)*(currentH/2) <= tensileYieldStressFrame ;

                    structurallyCompliant = all([directStressCompliant, shearStressCompliant, bendingStressCompliant]);

                    if structurallyCompliant
                        success = [success; i, j, k, l, currentFrameArea];
                    end
             
                end
            end
        end
    end


    successfulH = hRange(success(:,1));
    successfulB = bRange(success(:,2));
    successfulC = cRange(success(:,3));
    successfulT = tRange(success(:,4));
    successfulArea = success(:,5);

    objectiveValue = successfulArea;
    [~, idxMin] = min(objectiveValue);

  
    figure;
    scatter3(successfulH, successfulB, successfulC, 50, objectiveValue, 'o', ...
         'MarkerFaceColor', 'flat', 'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);

    xlabel('Web Height (H)');
    ylabel('Bottom Flange Width (B)');
    zlabel('Top Flange Width (C)');
    colorbar;
    title('Heavy Frame Omega-Section Optimisation');

    hold on;
    plot3(successfulH(idxMin), successfulB(idxMin), successfulC(idxMin),...
        'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    legend('Design Points', 'Minimum Objective');


    h = successfulH(idxMin);
    b = successfulB(idxMin);
    c = successfulC(idxMin);
    t = successfulT(idxMin);

   % plotOmegaFrame(h, b, c, t);
end
