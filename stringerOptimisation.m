function [skinThickness,stringerThickness,h,L,stringerSpacing] = stringerOptimisation(density)
      
    %stringer
            tensileYieldStressStringer = 490*1e+06;
            shearYieldStressStringer = 310*1e6;
            EStringer = 71*1e9;  

            density = 2780;
            L_eff = 0.5;         
            Ks = 5;       
            nu = 0.33;        
            b=0.02;

            %skin
            tensileYieldStressSkin = 331*1e6;
            shearYieldStressSkin = 290*1e6;
            ESkin = 73e9;      
         

            E = 71e9;
            shearYieldStress = 310*1e6;
            tensileYieldStress = 490*1e+06;

% skin thickness, spacing,str thick,h,l

    lb = [0.004, 0.1, 0.001, 0.01, 0.01];
    ub = [0.005, 0.5, 0.003, 0.045, 0.03];
    
    skinThicknessRange = lb(1):0.00005:ub(1);
    stringerSpacingRange = lb(2):0.005:ub(2);
    stringerThicknessRange = lb(3):0.00005:ub(3);
    hRange = lb(4):0.0005:ub(4);
    LRange = lb(5):0.0005:ub(5);

    [~, ~, ~, ~, totalSFlow] = plotShearFlow(shearYieldStressSkin);
    
    fprintf('Total full iterations to try: %d\n', length(skinThicknessRange));
    success = [];
    for i = 1:length(skinThicknessRange)
        disp(['Starting full iteration index: ', num2str(i)]);
        currentSkinThickness = skinThicknessRange(i);
        for j = 1:length(stringerSpacingRange)
            currentStringerSpacing = stringerSpacingRange(j);
            [~, y, numStringers] = stringerPlot(currentStringerSpacing);
    
            for k = 1:length(stringerThicknessRange)
                % Use the current stringer thickness from the range
                currentStringerThickness = stringerThicknessRange(k);
    
                for l = 1:length(hRange)
                    currentH = hRange(l);
    
                    for m = 1:length(LRange)
                        currentL = LRange(m);
    
                        currentStringerArea = ZStringerArea(currentStringerThickness, currentH, currentL);                   
                        directStress = stressPlot(currentStringerThickness, currentStringerArea, y, numStringers, density);
                        %[stringerStressCompliant, stringerEulerCompliant] = checkStringer(...
                           % directStress, tensileYieldStressStringer, EStringer, b, currentH, currentStringerArea, L_eff);
                        [stringerStressCompliant,stringerEulerCompliant] = checkStringer(directStress, tensileYieldStressStringer, ...
                            EStringer, currentH, currentStringerArea, currentL, currentStringerThickness);


                        [skinYieldCompliant, skinBucklingCompliant] = checkSkinBay(...
                            totalSFlow, currentSkinThickness, currentStringerSpacing, shearYieldStressSkin, Ks, ESkin, nu);
    
                        structurallyCompliant = all([stringerStressCompliant, stringerEulerCompliant, skinYieldCompliant, skinBucklingCompliant]);
    
                        if structurallyCompliant
                            success = [success; i, j, k, l, m, currentStringerArea, numStringers];
    
                        end
    
                    end 
                end 
            end 
        end 
        
    end
    
    
    successfulSkinThickness   = skinThicknessRange(success(:,1));
    successfulStringerSpacing = stringerSpacingRange(success(:,2));
    successfulStringerThickness = stringerThicknessRange(success(:,3));
    successfulH = hRange(success(:,4));
    successfulL = LRange(success(:,5)); 
    
    objectiveValue = success(:,6).*success(:,7);
    
    figure;
    scatter3(successfulSkinThickness, successfulStringerThickness, successfulH, 50, objectiveValue, 'filled');
    xlabel('Skin Thickness');
    ylabel('Stringer Thickness');
    zlabel('H');
    colorbar;
    title('Design Space: Stringer Area Objective');

    [minValue, idxMin] = min(objectiveValue);
    hold on;
    plot3(successfulSkinThickness(idxMin), successfulStringerThickness(idxMin), successfulH(idxMin),...
      'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    legend('Design Points', 'Minimum Objective');

    skinThickness = successfulSkinThickness(idxMin);
    stringerThickness = successfulStringerThickness(idxMin);
    h = successfulH(idxMin);
    t_s = successfulStringerThickness(idxMin);
    L = successfulL(idxMin);
    stringerSpacing = successfulStringerSpacing(idxMin);

end