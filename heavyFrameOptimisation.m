function [h,l,t] = heavyFrameOptimisation(maxN,maxS,maxM)
     
    tensileYieldStressFrame = 375*1e6;
    shearYieldStressFrame = 200*1e6;

    % lb = [ 0.02, 0.01,0.001,0.002 ];
    % ub = [ 0.06, 0.04, 0.004, 0.004 ];

    lb = [ 0.01, 0.01,0.01];
    ub = [ 0.5, 0.1, 0.1];
    
    hRange = lb(1):0.005:ub(1);
    lRange = lb(2):0.003:ub(2);
    tRange = lb(3):0.003:ub(3);


    success = [];
    for i = 1:length(hRange)
        currentH = hRange(i);

        for j = 1:length(lRange)
            currentL = lRange(j);

            for k = 1:length(tRange)
                currentT = tRange(k);
 
                    [currentFrameArea,currentFrameIxx] = IFrameProps(currentT, currentH, currentL); 

                    if maxN <= tensileYieldStressFrame*currentFrameArea
                        directStressCompliant = true;
                    else
                        directStressCompliant = false;
                    end
                         %apparently the web takes most the shear flow
                    if maxS <= shearYieldStressFrame * currentH*currentT
                        shearStressCompliant = true;
                    else
                        shearStressCompliant = false;
                    end
                    
                    if maxM <= tensileYieldStressFrame*currentFrameIxx/(currentH/2)
                        bendingStressCompliant = true;
                    else
                        bendingStressCompliant = false;
                    end

                    structurallyCompliant = all([directStressCompliant, shearStressCompliant, bendingStressCompliant]);
           
                    if structurallyCompliant
                            success = [success; i, j, k, currentFrameArea];
    
                    end
            
            end 
        end 
     %   disp(['Completed index: ' num2str(i)]);
    end

     
    successfulH   = hRange(success(:,1));
    successfulL = lRange(success(:,2));
    successfulT = tRange(success(:,3));
    successfulArea = (success(:,4));

    objectiveValue = successfulArea;
    [~, idxMin] = min(objectiveValue);
    % 
    % figure;
    % scatter3(successfulH, successfulL, successfulT, 50, objectiveValue, 'filled');
    % xlabel('H');
    % ylabel('L');
    % zlabel('T');
    % colorbar;
    % title('Heavy Frame Area Optimisation');
    % 
    % hold on;
    % plot3(successfulH(idxMin), successfulL(idxMin), successfulT(idxMin),...
    %   'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % legend('Design Points', 'Minimum Objective');

    h = successfulH(idxMin);
    l = successfulL(idxMin);
    t = successfulT(idxMin);



    
end