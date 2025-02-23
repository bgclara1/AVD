function [yieldCompliant, bucklingCompliant] = checkSkinBay(totalSFlow, skinThickness, StringerSpacing, shearYieldStress, Ks, E, nu)
   
    maxShearFlow = max(abs(totalSFlow));
    shearStress_skin = maxShearFlow / skinThickness;
    yieldCompliant = shearStress_skin <= shearYieldStress;
    tau_cr = Ks * (pi^2 * E) / (12 * (1 - nu^2)) * (skinThickness / StringerSpacing)^2;
    bucklingCompliant = shearStress_skin <= tau_cr;
    
end
