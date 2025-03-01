function [hoopStress, longitudinalStress, sphericalStress] = pressureStresses(D, skinThickness,pressureLoad,atmosphericPressure)

    hoopStress = (D/(2*skinThickness))*pressureLoad*atmosphericPressure;
    longitudinalStress = (D/(4*skinThickness))*pressureLoad*atmosphericPressure;
    sphericalStress = (D/(4*skinThickness))*pressureLoad*atmosphericPressure;

end