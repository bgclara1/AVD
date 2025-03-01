% MATLAB script to read, plot SC-0712 airfoil, and calculate centroid & flexural axis

clear
clc
close all


% Load the airfoil data
filename = 'sc20712.dat.txt';
data = load(filename);

% Extract x and y coordinates
x = data(:,1);
y = data(:,2);

% Compute the centroid (geometric center)
xc = mean(x);
yc = mean(y);

% % Plot the centroid
% plot(xc, yc, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');


% % Estimate the flexural axis (assuming it is close to 40% of the centroid location)
% x_flexural = 0.4 * xc;
% y_flexural = 0.4 * yc;
% plot(x_flexural, y_flexural, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% legend('Airfoil Shape', 'Centroid', 'Flexural Axis');
% 
% % Display centroid and flexural axis coordinates
% disp(['Centroid: (', num2str(xc), ', ', num2str(yc), ')']);
% disp(['Flexural Axis: (', num2str(x_flexural), ', ', num2str(y_flexural), ')']);

% Define the wingbox coordinates
x_wb = [x(23), x(83), x(124), x(184), x(23)]; % Front Spar, Rear Spar
y_wb = [y(23), y(83), y(124), y(184), y(23)]; % 0.2, 0.8

% x_wb = [x(18), x(83), x(124), x(189), x(18)]; % 0.15, 0.8
% y_wb = [y(18), y(83), y(124), y(189), y(18)];

% x_wb = [x(18), x(78), x(129), x(189), x(18)]; % 0.15, 0.75
% y_wb = [y(18), y(78), y(129), y(189), y(18)];


% Calculate wingbox centroid using the polygon method
A = 0;
Cx = 0;
Cy = 0;
for i = 1:length(x_wb)-1
    common_term = x_wb(i)*y_wb(i+1) - x_wb(i+1)*y_wb(i);
    A = A + common_term;
    Cx = Cx + (x_wb(i) + x_wb(i+1)) * common_term;
    Cy = Cy + (y_wb(i) + y_wb(i+1)) * common_term;
end
A = A / 2;
Cx = Cx / (6*A);
Cy = Cy / (6*A);

% Display the wingbox centroid
disp(['Wingbox Centroid: (', num2str(Cx), ', ', num2str(Cy), ')']);

ave_h = ((y(23) - y(184)) +  (y(83) - y(124)))/2

y_wb_ideal = [ave_h/2, ave_h/2, -ave_h/2, -ave_h/2, ave_h/2]; % 0.2, 0.8


% Plot the airfoil shape
figure;
hold on
plot(x, y, 'b-', 'LineWidth', 1.5);
plot(x_wb, y_wb, 'r--', 'LineWidth', 2);
plot(x_wb, y_wb_ideal, 'k-', 'LineWidth', 2);
plot(Cx, Cy, 'ko', 'MarkerSize', 8.5, 'MarkerFaceColor', '[0.6510    0.6510    0.6510]');
%plot(xc, yc, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend('Airfoil Shape', 'Real Wingbox', 'Idealised Wigbox', 'Shear Centre');
grid on;
axis equal;
xlabel('Chordwise Position (x)');
ylabel('Vertical Position (y)');
%title('SC-0712 Airfoil Shape');



% Calculate moment of inertia (Ix, Iy) using the parallel axis theorem
Ix = 0;
Iy = 0;
for i = 1:length(x_wb)-1
    common_term = x_wb(i)*y_wb(i+1) - x_wb(i+1)*y_wb(i);
    Ix = Ix + (y_wb(i)^2 + y_wb(i)*y_wb(i+1) + y_wb(i+1)^2) * common_term;
    Iy = Iy + (x_wb(i)^2 + x_wb(i)*x_wb(i+1) + x_wb(i+1)^2) * common_term;
end
Ix = Ix / 12;
Iy = Iy / 12;

% Display moment of inertia
disp(['Moment of Inertia (Ix): ', num2str(Ix)]);
disp(['Moment of Inertia (Iy): ', num2str(Iy)]);



% % Define the indices for the leading edge and wing box boundary
% start_idx = 23;
% end_idx = 184;
% 
% x2 = [x(start_idx:-1:1); x(end:-1:end_idx)];
% y2 = [y(start_idx:-1:1); y(end:-1:end_idx)];
% 
% % Compute the cumulative distance along the airfoil
% arc_length = 0;
% for i = 1:length(x2)-1
%     dx = x2(i+1) - x2(i);
%     dy = y2(i+1) - y2(i);
%     arc_length = arc_length + sqrt(dx^2 + dy^2);
% end
% 
% disp(['Distance along the airfoil from index ', num2str(start_idx), ' to ', num2str(end_idx), ': ', num2str(arc_length)]);
% 
% % Scaling for a linear chord distribution
% % Assume a new chord length C_new and original chord length C_old
% C_old = 1;  % Original chord length
% C_new = 3.0;              % Example: Normalize to unit chord
% scale_factor = C_new / C_old;
% 
% % Apply scaling to arc length
% scaled_arc_length = arc_length * scale_factor;
% disp(['Scaled arc length for linear chord distribution: ', num2str(scaled_arc_length)]);
% 
% 
% % Estimate Radius of Curvature of LE
% radius = 0.02;
% theta = linspace(0, 2*pi, 100);
% circle_x = (radius * cos(theta)) + radius;
% circle_y = radius * sin(theta);
% plot(circle_x, circle_y, 'k--', 'LineWidth', 1.5);
% 
% radius2 = 4;
% theta = linspace(0, 2*pi, 100);
% circle_x2 = (radius2 * cos(theta)) + 0.57;
% circle_y2 = radius2 * sin(theta) - radius2 + 0.07;
% plot(circle_x2, circle_y2, 'g-.', 'LineWidth', 1.5);
% xlim([-0.5,1.5])
% ylim([-0.2,0.2])
