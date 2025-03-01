clear
clc
close all


%% === Planform === %


% Given data
b = 65;         % Span in meters
sweep_c4 = 31.5;

% Calculate reference area and chord lengths
c_root = 13.6;
c_tip = 2.16;

% Convert quarter-chord sweep angle to radians
sweep_c4_rad = deg2rad(sweep_c4);

% Calculate Leading Edge (LE) and Trailing Edge (TE) sweep angles in radians
sweep_LE_rad = atan(tan(sweep_c4_rad) + (c_root - c_tip) / b);
sweep_TE_rad = atan(tan(sweep_c4_rad) - (c_root - c_tip) / b);

% Calculate x-positions of leading and trailing edges at root and tip
x_le_root = 0;                             % Leading edge at root (0, 0)
x_le_tip = (b / 2) * tan(sweep_LE_rad);    % Leading edge at tip
x_te_root = x_le_root + c_root;            % Trailing edge at root
x_te_tip = x_le_tip + c_tip;               % Trailing edge at tip

% Plot the wing planform
%figure;
hold on;
plot([0,b/2],[x_le_root,x_le_tip], 'k-', 'LineWidth', 2); % Leading edge
plot([0,b/2], [x_te_root,x_te_tip], 'k-', 'LineWidth', 2); % Trailing edge
plot([0,0], [x_le_root, x_te_root], 'k-', 'LineWidth', 2); % Root chord
plot([b/2,b/2], [x_le_tip, x_te_tip], 'k-', 'LineWidth', 2); % Tip chord
axis ij


% Add labels and title
xlabel('Spanwise Position [m]');
ylabel('Streamwise positing from wing root [m]');
%title('Wing Planform');
axis equal;
grid on;
xlim([0,35])
ylim([0,28])




%% === Fuselage === %%

fuselage_diam = 6.49;   % fuselage diameter (m)
fuselage_rad = fuselage_diam/2; % fuselage radius (m)

xline(fuselage_rad,'k--',LineWidth=2)



%% === SPARS === %%

% Calculate Front and Rear Spar sweep angles in radians
sweep_FrontSpar_rad = atan(tan(sweep_c4_rad) + 0.2*(c_root - c_tip) / b);
sweep_RearSpar_rad = atan(tan(sweep_c4_rad) - 0.8*(c_root - c_tip) / b);

% Calculate x-positions of leading and trailing edges at root and tip
x_fs_root = 0.2*c_root;                             % Front Spar at root
x_fs_tip = x_le_tip + 0.2 * c_tip;    % Front Spar at tip
x_rs_root = 0.8 * c_root;            % Rear Spar at root
x_rs_tip = x_le_tip + 0.8 * c_tip;               % Rear Spar at tip

plot([0,b/2],[x_fs_root,x_fs_tip], 'b-', 'LineWidth', 2.2); % Leading edge
plot([0,b/2], [x_rs_root,x_rs_tip], 'b-', 'LineWidth', 2.2); % Trailing edge





%% === Ribs === %%

%Rib_Locations = [4,5,6,7.2,8.5,10,12,15,19,24,30];
Rib_Locations = [4.81741567509198	6.38949949290096	7.96352337233209	9.54286442552014	11.1278270411260	12.7204513838344	14.3305523544488	15.9606484058625	17.6131554263673	19.2899345429857	20.9913206370639	22.7260840110450	24.5343144078981	26.4848207410255	28.2951722981210	30.1601667409953];


% Define two points (x1, y1) and (x2, y2)
x_fs = [0, b/2];  % X-coordinates of two points
y_fs = [x_fs_root, x_fs_tip];  % Y-coordinates of two points
% x_rs = [0, b/2];  % X-coordinates of two points
% y_rs = [x_rs_root, x_rs_tip];  % Y-coordinates of two points

% Interpolate y-values using linear interpolation
front_spar_interp = interp1(x_fs, y_fs, Rib_Locations, 'linear');
%rear_spar_interp = interp1(x_rs, y_rs, Rib_Locations, 'linear');

% Rib Line
%scatter(Rib_Locations,front_spar_interp)

% Define constants
m1 = 0.5072;      % Slope of first line
c1 = 10.88;       % Intercept of first line
m2 = 625/449;     % Slope of second line

% Calculate intersection points
x_rear_spar_rib = ((m2 .* Rib_Locations) + front_spar_interp - c1) ./ (m1 + m2);
y_rear_spar_rib = m1 * x_rear_spar_rib + c1;



% Plot Ribs

for i = 1:length(Rib_Locations)
    plot([Rib_Locations(i),x_rear_spar_rib(i)],[front_spar_interp(i),y_rear_spar_rib(i)],'r-',  'LineWidth', 2)
end



%% === Pseudo Ribs === %%

Pseudo_rib_Locations = linspace(fuselage_rad,b/2,25);

% Define two points (x1, y1) and (x2, y2)
x_le = [0, b/2];  % X-coordinates of two points
y_le = [x_le_root, x_le_tip];  % Y-coordinates of two points
% x_rs = [0, b/2];  % X-coordinates of two points
% y_rs = [x_rs_root, x_rs_tip];  % Y-coordinates of two points

% Interpolate y-values using linear interpolation
le_interp = interp1(x_le, y_le, Pseudo_rib_Locations, 'linear');
%rear_spar_interp = interp1(x_rs, y_rs, Rib_Locations, 'linear');

% Rib Line
%scatter(Rib_Locations,front_spar_interp)

% Define constants
m1 = 0.7184;      % Slope of first line
c1 = 2.72;       % Intercept of first line
m2 = 625/449;     % Slope of second line

% Calculate intersection points
x_pseudo_rib = ((m2 .* Pseudo_rib_Locations) + le_interp - c1) ./ (m1 + m2);
y_pseudo_rib = m1 * x_pseudo_rib + c1;


% Plot Pseudo Ribs

for i = 2:length(Pseudo_rib_Locations)
    plot([Pseudo_rib_Locations(i),x_pseudo_rib(i)],[le_interp(i),y_pseudo_rib(i)],'Color','[0.3922    0.8314    0.0745]',  'LineWidth', 2)
end




%% === Stringers === %%

taper = c_tip/c_root;
c_fuselage = c_root * (1 - (1-taper)*(2*fuselage_rad/b));

x_stringer_start = 0.7184*fuselage_rad + 2.72;
x_stringer_end = 0.5072*fuselage_rad + 10.88;
Stringer_Locations = linspace(x_fs_root,x_rs_root,20);
%Stringer_Locations = linspace(x_stringer_start,x_stringer_end,20);


% Stringer Line
%scatter(zeros(size(Stringer_Locations)),Stringer_Locations)

% Define constants
m1 = 0.7184;      % Slope of first line
c1 = 2.72;       % Intercept of first line
m2 = 0.5072;     % Slope of second line

% Calculate intersection points
x_stringer = (Stringer_Locations - 2.72)./(0.7184-0.5072);
y_stringer = m1 * x_stringer + c1;
x_stringer2 = x_stringer;
y_stringer2 = y_stringer;
Stringer_Locations2 = Stringer_Locations;


% Plot Stringers

for i = 1:length(Stringer_Locations)-1
    if x_stringer(i) > b/2
        x_stringer(i) = b/2;
        y_stringer(i) = 0.5072 * x_stringer(i) + Stringer_Locations(i);
    % elseif x_stringer(i) < fuselage_rad
    %     x_stringer(i) = fuselage_rad;
    %     y_stringer(i) = 0.5072 * x_stringer(i) + Stringer_Locations(i);
    end

    if x_stringer(i) > fuselage_rad
        x_stringer2(i) = fuselage_rad;
        y_stringer2(i) = 0.5072 * x_stringer2(i) + Stringer_Locations(i);
        Stringer_Locations2(i) = Stringer_Locations(i);
    end

    plot([0,x_stringer(i)],[Stringer_Locations(i),y_stringer(i)],'Color','[0.6510    0.6510    0.6510]',  'LineWidth', 1)
    plot([0,x_stringer2(i)],[Stringer_Locations2(i),y_stringer2(i)],'Color','[1,1,1]',  'LineWidth', 1)
   
end
%scatter(x_stringer,y_stringer)
%%
legend('','','','','Fuselage','Spars','', 'Ribs','','','','','','','','','','','','','', '','','','','','','','','','','', 'Pseudo Ribs','','','','','','','','','','','','','','','','','','','','','','', 'Stringers')




