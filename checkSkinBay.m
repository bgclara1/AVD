function [yieldCompliant, bucklingCompliant] = checkSkinBay(totalSFlow, skinThickness, StringerSpacing, shearYieldStress, Ks, E, nu)
   
    %[skinThickness, StringerSpacing]= [0.01,0.2]

% increase skin thickness to decrease shearStress_skin
% decrease stringer spaceing increases tau crit

    maxShearFlow = max(abs(totalSFlow));
    shearStress_skin = maxShearFlow / skinThickness;
    yieldCompliant = shearStress_skin <= shearYieldStress;
    
    tau_cr = Ks * (pi^2 * E) / (12 * (1 - nu^2)) * (skinThickness / StringerSpacing)^2;
    print = 'in function';
    bucklingCompliant = shearStress_skin <= tau_cr;


end

