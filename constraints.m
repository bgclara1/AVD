function [c, ceq] = constraints(x, shearYieldStress, tensileYieldStress, E, L_eff, Ks, nu, r)
    skinThickness = x(1);
    StringerSpacing = x(2);
    stringerThickness = x(3);
    h = x(4);
    L = x(5);

    numStringers = floor((2*pi*r)/StringerSpacing);
    stringerArea = ZStringerArea(stringerThickness, h, L);
    
    % Recalculate stresses and loads based on variables
    [~, ~, ~, ~, totalSFlow] = plotShearFlow(shearYieldStress);
    [~, y, ~] = stringerPlot(StringerSpacing);
    directStress = stressPlot(skinThickness, stringerArea, y, numStringers, 2780);

    % Stringer direct stress constraint
    stringerYieldCheck = max(abs(directStress)) - tensileYieldStress;

    % Euler buckling constraint (P_cr = pi^2*E*I/(L_eff^2))
    I = (stringerThickness * h^3)/12; % Simplified stringer I
    P_cr = (pi^2 * E * I)/(L_eff^2);
    sigma_cr_euler = P_cr/stringerArea;
    stringerEulerCheck = max(abs(directStress)) - sigma_cr_euler;

    % Skin bay yield constraint
    maxShearFlow = max(abs(totalSFlow));
    shearStress_skin = maxShearFlow / skinThickness;
    skinYieldCheck = shearStress_skin - shearYieldStress;

    % Skin bay buckling constraint (local plate buckling)
    tau_cr = Ks * (pi^2 * E) / (12 * (1 - nu^2)) * (skinThickness / StringerSpacing)^2;
    skinBucklingCheck = shearStress_skin - tau_cr;

    % Combine all constraints into a vector (c ≤ 0)
    c = [stringerYieldCheck;
         stringerEulerCheck;
         skinYieldCheck;
         skinBucklingCheck];

    ceq = []; % no equality constraints
end
