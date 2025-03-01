function [c, ceq] = constraints(x, shearYieldStressSkin, tensileYieldStressStringer, ESkin,EStringer, L_eff, Ks, nu, r)


    skinThickness = x(1);
    StringerSpacing = x(2);
    stringerThickness = x(3);
    h = x(4);
    L = x(5);

    stringerArea = ZStringerArea(stringerThickness, h, L);


    [~, ~, ~, ~, totalSFlow] = plotShearFlow(shearYieldStressSkin);
    [~, y, numStringers] = stringerPlot(StringerSpacing);
    directStress = stressPlot(skinThickness, stringerArea, y, numStringers);

    [stringerStressOK, stringerEulerOK] = checkStringer(directStress, tensileYieldStressStringer, EStringer, stringerThickness, h, stringerArea, L_eff);
    [skinYieldOK, skinBucklingOK] = checkSkinBay(totalSFlow, skinThickness, StringerSpacing, shearYieldStressSkin, Ks, ESkin, nu);

   
    c(1) = -double(stringerStressOK) + 0.5; 
    c(2) = -double(stringerEulerOK) + 0.5;   
    c(3) = -double(all(skinYieldOK)) + 0.5;     
    c(4) = -double(all(skinBucklingOK)) + 0.5;  

    ceq = [];
end
