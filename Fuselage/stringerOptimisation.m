function [skinThickness, stringerThickness, h, L, stringerSpacing] = optimizeStringerDesign(density)


    set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',25)
set(groot,'defaulttextfontsize',25)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')

    tensileYieldStressStringer = 490e6;  
    shearYieldStressStringer = 310e6;     
    EStringer = 71e9;                     

    tensileYieldStressSkin = 331e6;       
    shearYieldStressSkin = 290e6;         
    ESkin = 73e9;                         

    L_eff = 0.5;          
    Ks = 5;       
    nu = 0.33;        
    b = 0.02;      

    %   t_skin,  bs ,t_string, h ,  l
    lb = [0.004, 0.1, 0.0016, 0.01, 0.01];  
    ub = [0.007, 0.4, 0.005, 0.08, 0.06]; 

    skinThicknessRange = linspace(lb(1), ub(1), 20);
    stringerSpacingRange = linspace(lb(2), ub(2), 35);
    stringerThicknessRange = linspace(lb(3), ub(3), 15);
    hRange = linspace(lb(4), ub(4), 25);
    LRange = linspace(lb(5), ub(5), 10);

    [~, ~, ~, ~, totalSFlow] = plotShearFlow(shearYieldStressSkin);

    success = [];

  fprintf('Total iterations: %d\n', length(skinThicknessRange));

    for i = 1:length(skinThicknessRange)
        disp(i)
        for j = 1:length(stringerSpacingRange)
            for k = 1:length(stringerThicknessRange)
                for l = 1:length(hRange)
                    for m = 1:length(LRange)

                        currentSkinThickness = skinThicknessRange(i);
                        currentStringerSpacing = stringerSpacingRange(j);
                        currentStringerThickness = stringerThicknessRange(k);
                        currentH = hRange(l);
                        currentL = LRange(m);

                        if (currentH / currentL < 1.5) || (currentH / currentL > 3)
                            continue; 
                        end
                        if (currentH / currentStringerThickness > 50)
                            continue; 
                        end
                        if (currentL / currentStringerThickness > 30)
                            continue;
                        end

                        [~, y, numStringers] = stringerPlot(currentStringerSpacing);
                        currentStringerArea = ZStringerArea(currentStringerThickness, currentH, currentL);

                        directStress = stressPlot(currentStringerThickness, currentStringerArea, y, numStringers, density);

                        [stringerStressCompliant, stringerEulerCompliant] = checkStringer(directStress, tensileYieldStressStringer, ...
                            EStringer, currentH, currentStringerArea, currentL, currentStringerThickness);

                        [skinYieldCompliant, skinBucklingCompliant] = checkSkinBay(totalSFlow, currentSkinThickness, ...
                            currentStringerSpacing, shearYieldStressSkin, Ks, ESkin, nu);

                        structurallyCompliant = all([stringerStressCompliant, stringerEulerCompliant, skinYieldCompliant, skinBucklingCompliant]);

                        if structurallyCompliant
                            success = [success; currentSkinThickness, currentStringerSpacing, currentStringerThickness, ...
                                                currentH, currentL, currentStringerArea, numStringers];
                        end
                    end
                end
            end
        end
    end

    successfulSkinThickness = success(:, 1);
    successfulStringerSpacing = success(:, 2);
    successfulStringerThickness = success(:, 3);
    successfulH = success(:, 4);
    successfulL = success(:, 5);
    successfulStringerArea = success(:, 6);
    successfulNumStringers = success(:, 7);

    skinMass = pi * (3.2425^2 - (3.2425 - successfulSkinThickness).^2) * density; 
    stringerMass = (successfulStringerArea .* successfulNumStringers) * density;
    objectiveValue = skinMass + stringerMass + successfulNumStringers; 

    figure;
    scatter3(successfulSkinThickness, successfulStringerThickness, successfulH, 50, objectiveValue, 'filled');
    xlabel('Skin Thickness (m)');
    ylabel('Stringer Thickness (m)');
    zlabel('Stringer Height (m)');
    colorbar;
    title('Successful Stringer Optimisation Combinations');

    [~, idxMin] = min(objectiveValue);
    hold on;
    plot3(successfulSkinThickness(idxMin), successfulStringerThickness(idxMin), successfulH(idxMin),...
          'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r'); 
    legend('Viable Points', 'Minimum Mass Point');

    skinThickness = successfulSkinThickness(idxMin)
    stringerThickness = successfulStringerThickness(idxMin)
    h = successfulH(idxMin)
    L = successfulL(idxMin)
    stringerSpacing = successfulStringerSpacing(idxMin)
    numStringers = successfulNumStringers(idxMin)

end
