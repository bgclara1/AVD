function mass = objective(x, density, r)
    skinThickness = x(1);
    spacing = x(2);
    t_s = x(3);
    h = x(4);
    L = x(5);
    
    C = 2*pi*r;
    numStringers = floor(C/spacing);
    stringerArea = ZStringerArea(t_s, h, L);
    totalArea = numStringers * stringerArea + C * skinThickness;
    mass = density * totalArea;
    
end

