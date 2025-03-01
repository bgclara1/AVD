function [hoopStress, longitudinalStress, sphericalStress, FatigueCycles] = pressureStresses(D, skinThickness, pressureLoad)


    hoopStress = (pressureLoad * D) / (2 * skinThickness); 
    longitudinalStress = (pressureLoad * D) / (4 * skinThickness); 
    sphericalStress = (pressureLoad * D) / (4 * skinThickness); 

    sigma_fatigue = 138e6;
    m_fatigue = 5; 
    
    FatigueCycles = ((sigma_fatigue / hoopStress) ^ m_fatigue);


end



