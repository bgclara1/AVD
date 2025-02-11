n = 3.75;
W = 330800*9.81;
Q = n*W;    %radial load
r = 6.39; %check
P = 0; %tangential load
T = 0; %torque

q = (Q*sin(phi))/(pi*r);
num = 190;
for i = 1:36
    phi(i) = num-10*i;
    if phi(i-1) == 0
        num = 360;
    end
end