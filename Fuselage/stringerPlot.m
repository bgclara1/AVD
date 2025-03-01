function [x,y,numStringers] = stringerPlot(StringerSpacing)

    n = 3.75;
    W = 330800 * 9.81;
    r = 3.195;    

    C = 2 * pi * r;
   % StringerSpacing = 0.5; % Change
    numStringers = floor(C / StringerSpacing);

    angle = linspace(0, 360, numStringers + 1);
    
    x = r * cosd(angle); 
    y = r * sind(angle); 
    
end