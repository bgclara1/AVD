close all;

n = 3.75;
W = 330800 * 9.81;
Q = n * W;    % Radial load
r = 3.195;     % Fuselage radius
P = 0;        % Tangential load
T = 0;        % Torque


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SHEAR FLOW + SKIN THICKNESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(1) = 180;
for i = 2:37
    phi(i) = phi(i-1) - 10;
    if phi(i-1) == 0
        phi(i) = 350;
    end
end


q = (Q * sin(deg2rad(phi))) / (pi * r); 
p = (P * cos(deg2rad(phi))) / (pi * r); 
t = (T + P * r) / (2 * pi * r^2); 

totalSFlow = q + p + t;

theta = deg2rad(phi); 

circle = ones(length(phi))*r*1e6;
PlotSFlow = [];
for i = 1:length(circle)
    if totalSFlow(i)+circle(i)<circle(i)
        PlotSFlow(i) = abs(totalSFlow(i)) + circle(i);
    else
        PlotSFlow(i) = (totalSFlow(i)) + circle(i);
    end 
end

%{
figure;
polarplot(theta,PlotSFlow,'b-','LineWidth', 3)
hold on;
polarplot(theta,circle,'r-', 'LineWidth', 2)
ax = gca; 
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
title('Shear Flow Distribution Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('Shear Flow', 'Fuselage', 'Location', 'best');
%}

maxTotSF = max(abs(totalSFlow));
shearYieldStress = 187.06*1e+06; %    CHANGE
skinThickness = maxTotSF/shearYieldStress;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MAX STRESS + FUSELAGE MASS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given Data
M = 1e+07;
D = 2 * r;
C = 2 * pi * r;
StringerSpacing = 0.5; % Change
numStringers = floor(C / StringerSpacing);
StringerArea = 0.005;
equivBoomArea = 15 * skinThickness^2;
SingleBoomArea = StringerArea + equivBoomArea;

angle = linspace(0, 360, numStringers + 1);

x = r * cosd(angle); 
y = r * sind(angle); 

areaxarm = SingleBoomArea .* y;
Iarray = SingleBoomArea .* y.^2;


figure;
scatter(x, y, 'filled');
grid on;
axis equal; 
margin = 0.3 * max(abs([x, y])); 
xlim([-max(abs(x)) - margin, max(abs(x)) + margin]);
ylim([-max(abs(y)) - margin, max(abs(y)) + margin]);
xlabel('X-axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y-axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Stringer Positions Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');



I = sum(Iarray);
c = max(abs(y));

sigma = (M*c)/(I*10^6);

A = StringerArea*numStringers+C*skinThickness;

density = 2.78e+03; %change
FuselageWeightPUL = density*A;

sigma = (M * c) ./ (I * 1e6);
directStress = sigma .* y / c; 

figure;
hold on;
grid on;
plot3(x, y, zeros(size(x)), 'r-', 'LineWidth', 2);
quiver3(x, y, zeros(size(x)), zeros(size(x)), zeros(size(x)), directStress, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
xlabel('z');
ylabel('y');
zlabel('Direct stress \sigma_z (MPa)');
title('Direct stress distribution around fuselage');
view(3);
legend('Fuselage cross-section', 'Direct stress at each stringer');
axis equal;

