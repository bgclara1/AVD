function [skinThicknessPressure,domeThickness] = pressureThicknesses(D,pressureLoad,atmosphericPressure,yieldStress, nu)
    skinThicknessPressure = D*atmosphericPressure*pressureLoad/(2*yieldStress);
    domeThicknessRatio = (2-nu)/(1-nu);
    domeThickness = skinThicknessPressure/domeThicknessRatio;
end