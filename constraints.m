function [c, ceq] = constraints(x, shearYieldStress, tensileYieldStress, E, L_eff, Ks, nu, r)

    % Unpack Variables
    skinThickness = x(1);
    StringerSpacing = x(2);
    stringerThickness = x(3);
    h = x(4);
    L = x(5);

    % Calculate Areas
    stringerArea = ZStringerArea(stringerThickness, h, L);

    % Get Stress Distributions
    [~, ~, ~, ~, totalSFlow] = plotShearFlow(shearYieldStress);
    [~, y, numStringers] = stringerPlot(StringerSpacing);
    directStress = stressPlot(skinThickness, stringerArea, y, numStringers);

    % Check constraints (single boolean per constraint)
    [stringerStressOK, stringerEulerOK] = checkStringer(directStress, tensileYieldStress, E, stringerThickness, h, stringerArea, L_eff);
    [skinYieldOK, skinBucklingOK] = checkSkinBay(totalSFlow, skinThickness, StringerSpacing, shearYieldStress, Ks, E, nu);

    % Constraints (negative means constraint satisfied)
    c(1) = -double(stringerStressOK) + 0.5;  % Must be true (≤0)
    c(2) = -double(stringerEulerOK) + 0.5;   % Must be true (≤0)
    c(3) = -double(all(skinYieldOK)) + 0.5;      % Ensure single scalar result
    c(4) = -double(all(skinBucklingOK)) + 0.5;   % Ensure single scalar result

    ceq = [];
end
