function mass = objective(x, density, r)
    skinThickness = x(1);
    StringerSpacing = x(2);
    stringerThickness = x(3);
    h = x(4);
    L = x(5);

    numStringers = floor((2*pi*r)/StringerSpacing);
    stringerArea = ZStringerArea(stringerThickness, h, L);
    totalArea = stringerArea*numStringers + (2*pi*r)*skinThickness;
    mass = totalArea * density; % kg/m
    
end

