% Housekeeping ============================================================
clear
clc
close all

% Formatting
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',24)
set(groot,'defaulttextfontsize',24)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')

% INERTIAl LOADING SECTION ================================================

% Aircraft variables
b = 65;     % wingspan (m)
semi_span = b/2;    % semi span (m)
fuselage_diam = 6.49;   % fuselage diameter (m)
fuselage_rad = fuselage_diam/2; % fuselage radius (m)
c_root = 13.6;  % root chord (m) 
c_tip = 2.16;   % tip chord (m)
fuel_start = fuselage_rad;     % start of fuel tanks (m)
fuel_end = 30;      % end of fuel tanks (m)
taper = c_tip/c_root;   % taper ratio
engine1 = 10.8;  % spanwise location of engine 1 (m)
engine2 = 20.5;  % spanwise location of engine 2 (m)
y_lg = 5.685;    % spanwise landing gear location (m) 11.37
MTOW = 330800;  % Max takeoff weight (kg)
MLW = 0.85*MTOW;    % Max landing weight (kg)

% Load factors and manoeuver speeds
n_ult = 3.75;   % ultimate load factor
n_ult_2 = -1.5;   % ultimate load factor

% Define component masses for the load distribution (all kg)
wing = 24237/2;
nacelle = 4618/4;
engine = 23239/4;
pneumatic = 259;
engine_controls = 1104;
main_lg = 5123;
fuel_systems = 936;
hydraulics = 271.9;
anti_icing = 757.9;
fuel = (176.95/2) * 0.804 * 1000;  %156760/2;

% Set a suitable discretization (m)
dy = 0.1;

% Create vector of spanwise points
y = 0:dy:semi_span;
y_fuel = fuel_start:dy:fuel_end;

% % Plot the chord distribution of the wing
% chord = trapezoidal(0,c_root,semi_span,c_tip,y); 
% figure()
% plot(y,chord,'k', 'DisplayName', 'Chord Distribution')
% hold on
% xline(fuselage_rad,'r--',LineWidth=2,Label='Fuselage',LabelVerticalAlignment='middle')
% xline(fuel_start,'b--',LineWidth=2,Label='Fuel Start',LabelVerticalAlignment='middle')
% xline(fuel_end,'b--',LineWidth=2,Label='Fuel End',LabelVerticalAlignment='middle')
% hold off
% xlabel('Spanwise location y (m)')
% ylabel('Chord (m)')

% Match the ratios of the wing structure to that of the chord (i.e. wing weight will be more where chord is greater)
wing_initial = dist_start(wing*9.81*n_ult,semi_span,taper);
wing_final = taper*wing_initial;
wing_dist = trapezoidal(0,-wing_initial,semi_span,-wing_final,y);

wing_initial_2 = dist_start(wing*9.81*n_ult_2,semi_span,taper);
wing_final_2 = taper*wing_initial_2;
wing_dist_2 = trapezoidal(0,-wing_initial_2,semi_span,-wing_final_2,y);

% Carry out the same process for the fuel and various sub-systems
fuel_and_systems = fuel + pneumatic + engine_controls + fuel_systems + hydraulics + anti_icing;
fuel_initial = dist_start(fuel_and_systems*9.81*n_ult,fuel_end-fuel_start,taper);
fuel_final = taper*fuel_initial;
fuel_dist = trapezoidal(fuel_start,-fuel_initial,fuel_end,-fuel_final,y_fuel);

fuel_initial_2 = dist_start(fuel_and_systems*9.81*n_ult_2,fuel_end-fuel_start,taper);
fuel_final_2 = taper*fuel_initial_2;
fuel_dist_2 = trapezoidal(fuel_start,-fuel_initial_2,fuel_end,-fuel_final_2,y_fuel);

% Interpolate fuel_dist onto the y vector (so they can be plotted on the same graph)
fuel_dist_interp = interp1(y_fuel, fuel_dist, y, 'linear', 0);
fuel_dist_interp_2 = interp1(y_fuel, fuel_dist_2, y, 'linear', 0);

% Add the engine + nacelle weights to the overall distribution
total_engine = (engine + nacelle)*9.81*n_ult;
total_engine_2 = (engine + nacelle)*9.81*n_ult_2;

% Add point loads for engines and landing gear
engine1_index = round(engine1/dy) + 1;
engine2_index = round(engine2/dy) + 1;
lg_index = round(y_lg/dy) + 1;

engine_load = total_engine / dy; % Spread over one spanwise discretization
point_loads = zeros(size(y));

point_loads(engine1_index) = -engine_load;
point_loads(engine2_index) = -engine_load;
point_loads(lg_index) = -main_lg*9.81*n_ult / dy;

engine_load_2 = total_engine_2 / dy; % Spread over one spanwise discretization
point_loads_2 = zeros(size(y));

point_loads_2(engine1_index) = -engine_load_2;
point_loads_2(engine2_index) = -engine_load_2;
point_loads_2(lg_index) = -main_lg*9.81*n_ult_2 / dy;

% Add the distributions together
total_dist = wing_dist + fuel_dist_interp + point_loads;
total_dist_2 = wing_dist_2 + fuel_dist_interp_2 + point_loads_2;

% % Plot the Individual Weight distributions
% figure()
% plot(y, wing_dist, 'r--', 'DisplayName', 'Wing Load Distribution')
% hold on
% plot(y, fuel_dist_interp, 'b--', 'DisplayName', 'Fuel Load Distribution')
% plot(y, point_loads, 'g--', 'DisplayName', 'Point Loads')
% plot(y, total_dist, 'k', 'DisplayName', 'Total Load Distribution')
% xlabel('Spanwise location y (m)')
% ylabel('Load Distribution (N/m)')
% title('Individual Inertia Components')
% legend()
% hold off
% 
% % Plot just the combined distribution
% figure()
% plot(y, total_dist, 'k', 'DisplayName', 'Total Load Distribution')
% xlabel('Spanwise location y (m)')
% ylabel('Load Distribution (N/m)')
% title('Load Distribution at Ultimate Load Factor n_{ult} = 3.75')
% legend()
% hold off


%% =========================================================================


% AERODYNAMIC LOADING SECTION =============================================

% Calculate the elliptic lift distribution
L_index = round(fuselage_rad/dy) + 1;

L0 = (4*n_ult*MTOW*9.81) / (pi*b);
L = L0.*sqrt(1-((y-fuselage_rad)./((b/2 - fuselage_rad))).^2);
%L = L0.*sqrt(1-((y-fuselage_rad)./((b/2 - fuselage_rad)./cosd(28.5))).^2);
L = L(L_index:end);

L0_2 = (4*n_ult_2*MTOW*9.81) / (pi*b);
L_2 = L0_2.*sqrt(1-((y-fuselage_rad)./(b/2 - fuselage_rad)).^2);
L_2 = L_2(L_index:end);


% % Plot the aerodynamic loading
% figure()
% plot(y,L,'k',DisplayName='Lift Distribution')
% xlabel('Spanwise Location y (m)')
% ylabel('Load Distribution (N/m)')
% legend()


% % Plot Aerodynamic Loading and Inertial Loading
% figure()
% hold on
% plot(y, total_dist, 'k', 'DisplayName', 'Total Load Distribution')
% plot(y_wing,L,'b',DisplayName='Lift Distribution')
% xlabel('Spanwise Location y (m)')
% ylabel('Load Distribution (N/m)')
% legend()
% title('Symmetric Flight Inertial and Aerodynamic Loads')
% hold off
% 
% 
% % Plot Aerodynamic Loading and Inertial Loading Clipped
% figure()
% hold on
% plot(y, total_dist, 'k', 'DisplayName', 'Total Load Distribution (n = 3.75)')
% plot(y_wing,L,'b',DisplayName='Lift Distribution (n = 3.75)')
% plot(y, total_dist_2, 'r--', 'DisplayName', 'Total Load Distribution (n = -1.5)')
% plot(y_wing,L_2,'m--',DisplayName='Lift Distribution (n = -1.5)')
% xlabel('Spanwise Location y (m)')
% ylabel('Load Distribution (N/m)')
% legend()
% %title('Symmetric Flight Inertial and Aerodynamic Loads (Clipped)')
% ylim([-0.5e6 0.3e6])
% hold off



%% =========================================================================


% GROUND LOAD SECTION =====================================================

V_vert = 3.05;      % Limit descent velocity according to CS 25.473 (m/s)
gear_stroke = 0.217;    % Distance of main lg movement upon landing (m)
n_landing = (1/9.81)*(0.5*V_vert^2/gear_stroke);    % Landing load factor

landing_load = MLW*(9.81/2)*n_landing;  % Landing load supported by one wing (N)

% Create vector to plot the landing load
lg_point_load = zeros(size(y));
landing_index = round(y_lg/dy) + 1;
lg_point_load(landing_index) = landing_load / dy;

% Match the ratios of the wing structure to that of the chord (i.e. wing weight will be more where chord is greater)
wing_initial_land = dist_start(wing*9.81*n_landing,semi_span,taper);
wing_final_land = taper*wing_initial_land;
wing_dist_land = trapezoidal(0,-wing_initial_land,semi_span,-wing_final_land,y);

% Carry out the same process for the fuel and various sub-systems
fuel_and_systems_land = fuel*0.85 + pneumatic + engine_controls + fuel_systems + hydraulics + anti_icing;
fuel_initial_land = dist_start(fuel_and_systems_land*9.81*n_landing,fuel_end-fuel_start,taper);
fuel_final_land = taper*fuel_initial_land;
fuel_dist_land = trapezoidal(fuel_start,-fuel_initial_land,fuel_end,-fuel_final_land,y_fuel);

% Interpolate fuel_dist onto the y vector (so they can be plotted on the same graph)
fuel_dist_interp_land = interp1(y_fuel, fuel_dist_land, y, 'linear', 0);

% Add the engine + nacelle weights to the overall distribution
total_engine_land = (engine + nacelle)*9.81*n_landing;

engine_load_land = total_engine_land / dy; % Spread over one spanwise discretization
point_loads_land = zeros(size(y));

point_loads_land(engine1_index) = -engine_load_land;
point_loads_land(engine2_index) = -engine_load_land;

% Plot the landing loads with the inertial loads
total_dist_land = -wing_dist_land - fuel_dist_interp_land + point_loads_land + lg_point_load;

% figure()
% plot(y, total_dist_land, 'k', 'DisplayName', 'Total Load Distribution')
% xlabel('Spanwise location y (m)')
% ylabel('Load Distribution (N/m)')
% legend()
% title('Nose off Landing')
% hold off
% 
% % Plot the clipped graphs for improved clarity
% figure()
% plot(y, total_dist_land, 'k', 'DisplayName', 'Total Load Distribution')
% xlabel('Spanwise location y (m)')
% ylabel('Load Distribution (N/m)')
% legend()
% ylim([-0.5e6 0.5e6])
% title('Nose off Landing (Clipped)')
% hold off


%% =========================================================================

% Plot All Load Cases Aerodynamic Loading and Inertial Loading Clipped

figure()
hold on
plot(y, total_dist, 'b', 'DisplayName', 'Inertial Load (n = 3.75)')
plot(y(L_index:end),L,'r--',DisplayName='Lift Distribution (n = 3.75)')
% plot(y, total_dist_2, '-.','color','r', 'DisplayName', ['Inertial Load Distribution, Symmetric', newline, 'Flight (n = -1.5)'])
% plot(y(L_index:end),L_2,'-.','color','[1.0000 0.4667 0]',DisplayName='Lift Distribution, Symmetric Flight (n = -1.5)')
plot(y, total_dist_land, 'g:', 'DisplayName', 'Landing Inertial Load (n = 2.18)')
xline(fuselage_rad,'k-',LineWidth=2,Label='Fuselage',LabelVerticalAlignment='bottom')
xlabel('Spanwise Location y (m)')
ylabel('Load Distribution (N/m)')
hleg = legend('show');
hleg.String(end) = []; % delete the last legend entry of the very last plot
%title('Symmetric Flight Inertial and Aerodynamic Loads (Clipped)')
ylim([-0.5e6 0.3e6])
hold off



%% =========================================================================

% Combined Load and Lift Distributions


%dist1 = total_dist + [zeros(1,L_index), L(L_index+1:end)];
dist375 = total_dist(L_index:end) + L;
dist150 = total_dist_2(L_index:end) + L_2;
dist_land = total_dist_land(L_index:end);
y_wing = y(L_index:end);



% figure
% hold on
% plot(y_wing,dist375)
% plot(y_wing,dist150)
% plot(y_wing,dist_land)





%% =========================================================================

% Shear Force and Bending Moment

% Use numerical integration to calculate the shear force and bending moment distributions

% Reverse the y and total_dist vectors for integration starting at the free
% end (see diagram on pg. 23 of slides to see how ordinate system is set up)
y_reversed = flip(y_wing);
dist375_rev = flip(dist375);
dist150_rev = flip(dist150);
dist_land_rev = flip(dist_land);

% Integrate to find Shear Force
sf375_reversed = cumtrapz(y_reversed, dist375_rev);
sf150_reversed = cumtrapz(y_reversed, dist150_rev);
sf_land_reversed = cumtrapz(y_reversed, dist_land_rev);

% Integrate to find Bending Moment
bm375_reversed = cumtrapz(y_reversed, sf375_reversed);
bm150_reversed = cumtrapz(y_reversed, sf150_reversed);
bm_land_reversed = cumtrapz(y_reversed, sf_land_reversed);


% Reverse the results back
sf375 = flip(sf375_reversed);
bm375 = flip(bm375_reversed);

sf150 = flip(sf150_reversed);
bm150 = flip(bm150_reversed);

sf_land = flip(sf_land_reversed);
bm_land = flip(bm_land_reversed);


% 
% % Plot Shear Force Distribution
% figure()
% hold on
% plot(y_wing,-sf375,'color','[0.6 0.8 1]',DisplayName='Symmetric Flight (n = 3.75)')
% plot(y_wing,-sf150,'-.','color','[1.0000 0.4667 0]',DisplayName='Symmetric Flight (n = -1.5)')
% plot(y_wing,-sf_land, '--','color','[0.3922 0.8314 0.0745]', 'DisplayName', 'Landing (n = 2.18)')
% xlabel('Spanwise location y (m)')
% ylabel('Shear Force (N)')
% title('Shear Force Distribution')
% legend()
% grid on
% hold off
% 
% % Plot Bending Moment Distribution
% figure()
% hold on
% plot(y_wing,bm375,'color','[0.6 0.8 1]',DisplayName='Symmetric Flight (n = 3.75)')
% plot(y_wing,bm150,'-.','color','[1.0000 0.4667 0]',DisplayName='Symmetric Flight (n = -1.5)')
% plot(y_wing,bm_land, '--','color','[0.3922 0.8314 0.0745]', 'DisplayName', 'Landing (n = 2.18)')
% xlabel('Spanwise location y (m)')
% ylabel('Bending Moment (N/m)')
% title('Bending Moment Distribution')
% legend()
% grid on
% hold off



%% =========================================================================

% TORQUE =============

V_A = 92.6;
V_D = 192.4;

% Pitch Mom
AOA_dist = -0.6 + (3/(b/2 - fuselage_rad)).*y(1:(length(y)-L_index+1)) + 3.7; % Setting + Twist + AOA Cruise
[CM_alpha,CM0] = CM_interp(AOA_dist);

wing_weight = wing_dist + point_loads;
fuel_weight = fuel_dist_interp;
wing_weight_2 = wing_dist_2 + point_loads_2;
fuel_weight_2 = fuel_dist_interp_2;
wing_weight_land = -wing_dist_land + point_loads_land + lg_point_load;
fuel_weight_land = -fuel_dist_interp_land;


T375 = torque(L,wing_weight,fuel_weight,V_A,y_wing,CM0,CM_alpha,engine1_index,engine2_index,L_index,c_root,c_tip,dy,y,false);
T150 = torque(L_2,wing_weight_2,fuel_weight_2,V_A,y_wing,CM0,CM_alpha,engine1_index,engine2_index,L_index,c_root,c_tip,dy,y,false);
T_land = torque(L,wing_weight_land,fuel_weight_land,V_A,y_wing,CM0,CM_alpha,engine1_index,engine2_index,L_index,c_root,c_tip,dy,y,true);


% % Plot Torque Distribution
% figure()
% hold on
% plot(y_wing,T375,'color','[0.6 0.8 1]',DisplayName='Symmetric Flight (n = 3.75)')
% plot(y_wing,T150,'-.','color','[1.0000 0.4667 0]',DisplayName='Symmetric Flight (n = -1.5)')
% plot(y_wing,T_land, '--','color','[0.3922 0.8314 0.0745]', 'DisplayName', 'Landing (n = 2.18)')
% xlabel('Spanwise location y (m)')
% ylabel('Torque (Nm)')
% title('Torque Distribution')
% legend()
% grid on
% hold off


%% FINAL LIMITING LOAD CASES ===============================================

% Shear Forces (N)
S_pos = max_of_vectors(abs(sf375),abs(sf_land));
S_neg = -sf150;

% Bending Moment (N/m)
M_pos = max_of_vectors(bm375,bm_land);
M_neg = bm150;

% Torque (Nm)
T_pos = max_of_vectors(T150,T_land);
T_neg = T375;

% figure()
% plot(y_wing,S_pos,'b',DisplayName='S^{+}')
% hold on
% plot(y_wing,S_neg,'r',DisplayName='S^{-}')
% xlabel('Spanwise location y (m)')
% ylabel('Shear (N)')
% legend()
% grid on
% hold off
% 
% figure()
% plot(y_wing,M_pos,'b',DisplayName='M^{+}')
% hold on
% plot(y_wing,M_neg,'r',DisplayName='M^{-}')
% xlabel('Spanwise location y (m)')
% ylabel('Bending Moment (Nm)')
% legend()
% grid on
% hold off
% 
% figure()
% plot(y_wing,T_pos,'b',DisplayName='T^{+}')
% hold on
% plot(y_wing,T_neg,'r',DisplayName='T^{-}')
% xlabel('Spanwise location y (m)')
% ylabel('Torque (Nm)')
% legend()
% grid on
% hold off

%% OPTIMISATION RESULTS ====================================================

% Vectors of stringer design parameters to test
t_s_vec = (3:1:10)*0.001;
b1_vec = 0.25:0.05:0.85;
h_vec = (50:5:100)*0.001;

% Store the results of the combination of each of these design variables
results = OptFunction(t_s_vec, b1_vec, h_vec, M_pos, sf375, T375, c_root, c_tip, y_wing);

%% RETAIN ONLY ACCEPTABLE COMBINATIONS =====================================

% Define limit values
Limit_sigma_compressive = 345e6;
tau_tresca = 175e6;

% Retain only results meeting criteria
filtered_results = results([]);
filtered_index = 1;

for i = 1:length(results)
    if results(i).max_combination <= 0.99 && ...
       results(i).min_efficiency >= 0.7 && ...
       results(i).max_sigma_crit < Limit_sigma_compressive && ...
       results(i).max_tau < tau_tresca
   
        filtered_results(filtered_index) = results(i);
        filtered_index = filtered_index + 1;
    end
end

acc_t_s = [filtered_results.t_s];
acc_b1 = [filtered_results.b1];
acc_h = [filtered_results.h];

norms = zeros(length(filtered_results),1);

for i = 1:length(filtered_results)
norms(i) = norm(filtered_results(i).obj(1:end-1));
end

%% STEP THESE RESULTS FOR MANUFACTURABILITY ================================

manufacturable_results = DetailedOptFunction(acc_t_s, acc_b1, acc_h, M_pos, sf375, T375, c_root, c_tip, y_wing);

% Extract values into arrays for vectorized filtering
max_combination = [manufacturable_results.max_combination];
min_efficiency = [manufacturable_results.min_efficiency];
max_sigma_crit = [manufacturable_results.max_sigma_crit];
max_tau = [manufacturable_results.max_tau];

% Logical indexing for filtering
valid_indices = (max_combination <= 0.99) & ...
                (min_efficiency >= 0.7) & ...
                (max_sigma_crit < Limit_sigma_compressive) & ...
                (max_tau < tau_tresca);

% Apply filtering
filtered_results = manufacturable_results(valid_indices);

final_t_s = [filtered_results.t_s];
final_b1 = [filtered_results.b1];
final_h = [filtered_results.h];

% Find most optimal solutions
[minNorms, minIdexes] = mink(norms,1);

% Combinations of values that correspond to these minimums
Design = filtered_results(minIdexes);

%% PLOT THE MOST OPTIMAL SOLUTION ==========================================

% Stop the rib thickness being rounded to 0
filtered_results(1).t_rib(filtered_results(1).t_rib == 0) = 1e-3;

figure()
plot(y_wing,filtered_results(1).t_stringers)
hold off

figure()
plot(y_wing,filtered_results(1).t_rib)
hold off
%% RUN ON FINAL GEOMETRY ===================================================
    
% Final run for the optimised geometry with the varying stringer spacing

E = 72e9;               % Youngs modulus (Pa)
rho = 2710;              % Density of Al (kg/m^3) 
M = M_pos;              % Bending moment (N/m)
c = linspace(c_root,c_tip,length(y_wing));      % Wing chord (m)
b2 = 0.0758.*c;             % Wing box height (m)
width = 0.6.*c;             % Wing box width  (m)

% Skin thickness with no stringers
t_no_stringers = ((M.*width)./(3.62.*E.*b2)).^(1/3);   % Skin thickness with no stringers (m)
Area_no_stringers = width.*t_no_stringers;          % Panel area with no stringers (m^2)

% Skin thickness with striners
b1 = 0.5;                     % Stringer pitch   (m)
n = width./b1;              % No. stringers       (/)
d_to_h = 0.3;               % width to height ratio (/)
h = 90e-3;                     % Stringer height         (m) 
d = h*d_to_h;               % Stringer width          (m)
t_s = 9e-3;                    % Stringer thickness      (m)
A_s = t_s*h + 2*t_s*d;      % Z-stringer area formula (m^2)
t_eff = ((M.*b1.^2)./(3.62.*E.*width.*b2)).^(1/3);       % Effective thickness (m)
N = M./b2./width;             % Axial load per unit length (N/m)
sigma = N./(t_eff);  % Critical global buckling stress (N/m^2)
t_stringers = t_eff - A_s./(b1);  % Skin thickness with stringers (m)
Area_with_stringers = n.*b1.*t_eff;             % Panel area with stringers (m^2)

% NOTE SOME OF THE ELEMENTS OF t_stringers WILL BE NEGATIVE WHEN A_s/b1 >
% t_eff. AS SUCH ANY ELEMENTS WHERE element - A_s/b1 < A_s/b1 ARE REPLACED 
% WITH THE VALUE OF 2*A_s/b1, THIS STOPS THE SHEAR STRESS EXPLODING AS 
% t_stringers APPROACHES 0. 

% Compute As/b1
threshold = A_s / b1;
% Apply the replacement condition
t_stringers((t_stringers - threshold) < threshold) = threshold;

t_stringers_man = ceil(1000.*t_stringers);
t_stringers_man = t_stringers_man./1000;

% UNCOMMENT THESE PLOTS TO COMPARE AREA OF STRINGER VS NO STRINGER
% Plot to compare affect of stringer vs. no stringers
% figure()
% plot(y_wing,t_no_stringers)
% hold on
% plot(y_wing,t_stringers)
% hold off
% 
% % Plot of different sectional area with and without stringers
% figure()
% plot(y_wing,Area_no_stringers,DisplayName='No Stringers')
% hold on
% plot(y_wing,Area_with_stringers,DisplayName='Stringers')
% legend()
% hold off


% Calculate the necessary ratios to apply the Catchpole correction factor
thickness_ratio = t_s./(t_stringers);
area_ratio = A_s./(b1.*t_stringers);

% Use the digitised Catchpole diagram to calculate the correction factor
% for the given vector of ratios
[stress_ratio, new_thickness] = catchpole_interpolation(thickness_ratio, area_ratio);

sigma_crit = stress_ratio.*sigma;        % Local critical buckling stress (N/m^2)

[efficiency, ~] = Farrar_interpolation(thickness_ratio, area_ratio);

rib_spacing_skin = ((efficiency./sigma_crit).^2).*N.*E;      % Optimal rib spacing (m) 

m_skinstringer = sum(2700.*t_eff.*0.1.*width);

%% SPAR DESIGN =============================================================

S = sf375;  % Shear force distribution (N)
T = -T375;   % Torque distribution (Nm)

q_torque = T./(2.*width.*b2);          % Shear flow due to torue (N/m)
q_shear = -S./(2.*b2);                  % Shear flow due the shear force (N/m)

q_fs = abs(q_torque + q_shear);          % Shear flow in front spar (N/m)
q_rs = abs(q_shear - q_torque);          % Shear flow in rear spar (N/m)

% Max Ks when a = b, can discretise the a distribution to ensure this
% constraint is met
a = b2;    % Set the vertical web spacing equal to the wing box height for now (m)

% Now start a cumulative spacing distribution such that at the start of the wing a = b2
% and the next spacing will be equal to the value of a at the span wise location = a etc. 
rib_spacing = calculate_stepwise_rib_spacing(y_wing, b2);

% Now Ks will be constant along the span and be equal to 13.4
Ks = 13.4;
Ks_test = 8.1;


% Calculate the thickness of the front and rear spars
t_front = (abs((q_fs.*b2)./(Ks*E))).^(1/3);
t_rear = (abs((q_rs.*b2)./(Ks*E))).^(1/3);
t_front_test = (abs((q_fs.*b2)./(Ks_test*E))).^(1/3);
t_rear_test = (abs((q_rs.*b2)./(Ks_test*E))).^(1/3);

% Mass of the spars
m_unstiffened_sapr = sum(t_front_test.*0.1.*2700);
m_spar_front = sum(t_front.*0.1.*2700) + 5e-3*2700*62*0.1;    % 62 = no. stiffeners, and 0.1m = average stiffener height
m_spar_rear = sum(t_rear.*0.1.*2700) + 5e-3*2700*62*0.1;

t_front_man = ceil(1000.*t_front);
t_front_man = t_front_man./1000;
t_rear_man = ceil(1000.*t_rear);
t_rear_man = t_rear_man./1000;

% Calculate the shear stress in the spars
tau_front = q_fs./t_front;
tau_rear = q_rs./t_rear;

% figure()
% plot(y_wing,t_front)
% hold on
% plot(y_wing,t_rear)
% hold off
% 
% figure()
% plot(y_wing,tau_front,"DisplayName","Front")
% hold on
% plot(y_wing,tau_rear,"DisplayName","Rear")
% legend()
% hold off

tau_skin = q_torque ./ t_stringers;             % Shear stress in skin (Pa)

% Set last value of tau equal to that of the penultimate to aoid dividing by 0
tau_skin(end) = tau_skin(end - 1);

tau_tresca = 175e6;                             % TRESCA shear yield stress (Pa)

combination = (tau_skin./tau_tresca).^2 + sigma./sigma_crit;

% figure()
% plot(y_wing,combination)
% hold off

%% RIB DESIGN =============================================================

I_skin_panel = (width.*(t_eff).^3)./12 + width.*t_eff.*(b2./2).^2;          % Moment of inertia (m^4)

Crush_load = 0.5.*((M.^2).*rib_spacing_skin.*b2.*(t_eff.*width))./(E.*I_skin_panel.^2);       % Crush load (N)

% Without Cutouts
t_rib_nocutouts = ((Crush_load.*(b2).^2)./(8.1*E*width)).^(1/3);           % Rib thickness (m)
t_rib_nc_man = ceil(1000.*t_rib_nocutouts);
t_rib_nc_man = t_rib_nc_man./1000;

% With Cutouts
t_rib = ((Crush_load.*(b2).^2)./(8.1/0.47*E*width)).^(1/3);           % Rib thickness (m)
t_rib_man = ceil(1000.*t_rib);
t_rib_man = t_rib_man./1000;

% Calculate masses
m_rib = 7*5e-3*2700*5 + 5*4e-3*2700*3 + 2*3e-3*2700*2 + 2*2e-3*2700*1;  
m_rib_nc = 2*7e-3*2700*5 + 5*6e-3*2700*5 + 5*5e-3*2700*3 + 2*4e-3*2700*2 + 2*3e-3*2700*1;   

sigma_rib = Crush_load./(t_rib.*width);         % Crush stress (Pa)


%% CREATE PLOTS ============================================================

figure()
plot(y_wing,t_stringers,'r')
hold on
plot(y_wing,t_stringers_man,'b--')
hold off
ylabel("Skin Thickness (m)")
xlabel("Wing Span (m)")
legend("Optimal Distribtuion","Maunfacturable")

figure()
plot(y_wing,t_rib,'r')
hold on
plot(y_wing,t_rib_man,'b--')
ylabel("Rib Thickness (m)")
xlabel("Wing Span (m)")
legend("Optimal Distribtuion","Maunfacturable")

figure()
plot(y_wing,efficiency,'r')
hold off
ylabel("Farrar Efficiency")
xlabel("Wing Span (m)")
ylim([0.5 1])

figure()
plot(y_wing,rib_spacing_skin,'b')
xlabel("Wing Span (m)")
ylabel("Optimal Rib Spacing (m)")

figure()
plot(y_wing,t_front_man,'b')
hold on
plot(y_wing,t_rear_man,'r')
xlabel("Wing Span (m)")
ylabel("Spar Thickness (m)")

figure()
plot(y_wing,combination,'b')
hold on
xlabel("Wing Span (m)")
ylabel("Combined Loading Value")

% % Define vectors
% t_s_vec = (3:1:10) * 0.001;  % 8 values
% b1_vec = 0.25:0.05:0.85;     % 13 values
% h_vec = (50:5:100) * 0.001;  % 11 values
% 
% % Create all combinations using meshgrid
% [T, B, H] = meshgrid(t_s_vec, b1_vec, h_vec);
% 
% % Convert to column vectors for scatter3
% T = T(:);
% B = B(:);
% H = H(:);

% % Plot the 3D scatter plot
% figure;
% scatter3(T, B, H, 50, 'filled'); % 50 is marker size, 'filled' for solid dots
% hold on
% % Plot acceptable solutions (red)
% scatter3(acc_t_s, acc_b1, acc_h, 80, 'r', 'filled'); % Red dots, slightly larger
% xlabel('t_{s} (m)');
% ylabel('b_{1} (m)');
% zlabel('h (m)');
% legend('All combinations', 'Acceptable solutions');
% grid on;
% view(3); 
% 
% figure;
% scatter3(acc_t_s, acc_b1, acc_h, 50, 'filled');
% hold on
% scatter3(final_t_s, final_b1, final_h, 80,'r', 'filled'); % Red dots, slightly larger
% hold on
% scatter3(filtered_results(1).t_s, filtered_results(1).b1, filtered_results(1).h, 100, 'g', 'filled'); % Green pentagon marker
% xlabel('t_{s} (m)');
% ylabel('b_{1} (m)');
% zlabel('h (m)');
% legend('Acceptable solutions','Manufacturable solutions','Optimal Solution');
% grid on;
% view(3); 

figure()
plot(y_wing,sigma_crit)
hold on
yline(345e6,'r--')
hold off
xlabel("Wing Span (m)")
ylabel("Stress (MPa)")



%%

plot(y_wing,rib_spacing_skin)
int_rib = interp1(y_wing,rib_spacing_skin,fuselage_rad+1.57)
x = fuselage_rad;
i = 1;

while x < 32.5
    rib_pos(i) = interp1(y_wing,rib_spacing_skin,x);
    x = x + rib_pos(i);
    i = 1 + i;
    rib_position(i) = x;
    
end

rib_position

%% FUNCTIONS ===============================================================

function rib_spacing = calculate_stepwise_rib_spacing(span_locations, chord)
    % Function to calculate stepwise rib spacing along the span of a wing.
    % The spacing remains constant between ribs and matches the dimensions of span_locations.
    %
    % Inputs:
    %   span_locations - Vector of spanwise locations along the wing
    %   chord - Vector of chord lengths corresponding to the span locations
    %
    % Output:
    %   rib_spacing - Vector of stepwise rib spacings along the wing span

    current_position = span_locations(1);   % Start at the first span location
    rib_positions = current_position;       % Store rib positions
    rib_values = [];                        % Store corresponding rib spacing values

    % Calculate rib positions and spacing values
    while current_position < max(span_locations)
        current_chord = interp1(span_locations, chord, current_position, 'linear');
        rib_values = [rib_values, current_chord];
        current_position = current_position + current_chord;
        if current_position <= max(span_locations)
            rib_positions = [rib_positions, current_position];
        end
    end

    % Initialize the final rib_spacing vector
    rib_spacing = zeros(size(span_locations));

    % Assign stepwise values along span_locations
    for i = 1:length(rib_positions)-1
        % Find indices within current rib section
        indices = span_locations >= rib_positions(i) & span_locations < rib_positions(i+1);
        rib_spacing(indices) = rib_values(i);
    end

    % Assign the last value to the remaining points
    rib_spacing(span_locations >= rib_positions(end)) = rib_values(end);

end



% Trapezoidal distribution function
function y = trapezoidal(x1, y1, x2, y2, x)
    % Inputs:
    %   x1, y1 - Coordinates of the first point
    %   x2, y2 - Coordinates of the second point
    %   x      - Vector of x-values where y-values are to be calculated
    %
    % Output:
    %   y      - Vector of corresponding y-values for the given x-values
    m = (y2 - y1) / (x2 - x1);
    c = y1 - m * x1;
    y = m * x + c;
end

function start = dist_start(weight,length,taper)
% Inputs:
    %   weight - weight of component being distributed
    %   length - length of distribution
    %   taper  - ratio of start and end distribution values (equivalent to taper ratio)
    %
    % Output:
    %   start  - start value of distribution
    start = (2*weight) / (length*(1+taper)); % (from A = 1/2(a+b)h)

end



function [CM_alpha,CM0] = CM_interp(AOA_dist)

% Load the data, skipping metadata
filename = 'SC0712.csv';
opts = detectImportOptions(filename);
opts.DataLines = [12 inf]; % Start reading from line 11 where the table begins
data = readtable(filename, opts);

% Extract relevant columns
AOA = data.Alpha; % Angle of Attack
Cm = data.Cm;     % Moment Coefficient

% Perform interpolation
CM_alpha = interp1(AOA, Cm, AOA_dist, 'spline');
CM0 = interp1(data.Cl, data.Cm, 0, 'linear');

% % Plot the results
% figure;
% plot(AOA, Cm, 'o', AOA_dist, CM_alpha, '-');
% xlabel('Angle of Attack (deg)');
% ylabel('Moment Coefficient Cm');
% title('Interpolated Cm vs AOA');
% legend('Original Data', 'Interpolated Data');
% grid on;

end



function load = dist2load(dist,dy)

    load = zeros(size(dist));
    for i = 1:length(dist)
        if i == length(dist)
            load(i) = (dist(i) + 0) * dy / 2;
        else
            load(i) = (dist(i) + dist(i+1)) * dy / 2;
        end
    end

end


function T = torque(L,inertia,fuel,V,y_wing,CM0,CM_alpha,engine1_index,engine2_index,L_index,c_root,c_tip,dy,y,land)

    c_bar = 8.75;
    c = linspace(c_root,c_tip,length(y_wing));
    S_ref = 469.44;
    rho0 = 1.225;
    x_sc = 0.50;

    % Pitch Mom
    
    L_tot = trapz(y_wing,L);
    M0 = (CM_alpha + CM0)*S_ref*c_bar*0.5*rho0*(V/3.6)^2;
    

    M0_dist = dist2load(M0 .* (L / (2*L_tot*dy)),dy);

    
    % Lift Moment
    if land == false
        x_ac = (x_sc - 0.25) * c;
    else
        x_ac = 0;
    end
    LM = L .* x_ac;
   
    % Thrust Moment
    thrust = 290070;
    z_T = 2.3;
    thrust_loads = zeros(size(y));
    thrust_loads(engine1_index) = thrust;
    thrust_loads(engine2_index) = thrust;
    
    TM = thrust_loads(L_index:end) * z_T;
    
    
    %Inertia Moment
    x_wing_max_thick = 0.37;
    x_w = -(x_sc - x_wing_max_thick) * c;
    IM = inertia(L_index:end) .* x_w;


    x_fuel = 0.48;
    x_f = -(x_sc - x_fuel) * c;
    FM = fuel(L_index:end) .* x_f;

    
    
    % Total Moments
    M = LM + TM + IM + FM + M0_dist;
    
    % 
    % figure()
    % hold on
    % plot(y_wing,LM)
    % plot(y_wing,TM)
    % plot(y_wing,IM)
    % plot(y_wing,M0_dist)
    % legend('LM','TM','IM','M0')
    
    
    % Torque
    M_rev = flip(M);
    torque_rev = cumtrapz(flip(y_wing),M_rev);
    T = flip(torque_rev);

end

function maxVector = max_of_vectors(varargin)
    % Get the number of vectors
    numVectors = nargin;
    
    % Check if at least one vector is provided
    if numVectors < 1
        error('At least one input vector is required.');
    end
    
    % Initialize the maxVector with the values of the first vector
    maxVector = varargin{1};
    
    % Iterate through all the vectors and find the max at each index
    for i = 2:numVectors
        maxVector = max(maxVector, varargin{i});
    end
end



function arc = D_Cell_Arc(c)

    % Load the airfoil data
    filename = 'sc20712.dat.txt';
    data = load(filename);
    
    % Extract x and y coordinates
    x = data(:,1);
    y = data(:,2);

    % Define the indices for the leading edge and wing box boundary
    start_idx = 23;
    end_idx = 184;
    
    x2 = [x(start_idx:-1:1); x(end:-1:end_idx)];
    y2 = [y(start_idx:-1:1); y(end:-1:end_idx)];
    
    % Compute the cumulative distance along the airfoil
    arc_length = 0;
    for i = 1:length(x2)-1
        dx = x2(i+1) - x2(i);
        dy = y2(i+1) - y2(i);
        arc_length = arc_length + sqrt(dx^2 + dy^2);
    end

    arc = arc_length .* c;

end



function Ks = curved_Ks(a,b_curve,r_DCell,t_DCell)
    % Determine Ks
    for i = 1:length(b_curve)
        if a > b_curve(i)
            a_b(i) = a./b_curve(i);
            b_rt(i) = b_curve(i) ./ sqrt(r_DCell(i).*t_DCell);
            %a_rt(i) = 0;
    
            if 8 < b_rt(i) && b_rt(i) < 9
    
                if  a_b(i) < 1.25
                    Ks(i) = 23.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 21.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 18.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 17;
                elseif 3 < a_b(i)
                    Ks(i) = 16;
                end
    
            elseif 7 < b_rt(i) && b_rt(i) < 8
    
                if  a_b(i) < 1.25
                    Ks(i) = 20.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 18;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 16.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 15;
                elseif 3 < a_b(i)
                    Ks(i) = 13.5;
                end
    
            elseif 6 < b_rt(i) && b_rt(i) < 7
    
                if  a_b(i) < 1.25
                    Ks(i) = 17;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 15.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 13.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 12.5;
                elseif 3 < a_b(i)
                    Ks(i) = 11;
                end
    
            elseif 5 < b_rt(i) && b_rt(i) < 6
    
                if  a_b(i) < 1.25
                    Ks(i) = 14.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 13;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 11.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 10.5;
                elseif 3 < a_b(i)
                    Ks(i) = 9.5;
                end
    
            elseif 3 < b_rt(i) && b_rt(i) < 5
    
                if  a_b(i) < 1.25
                    Ks(i) = 11;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 10;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 8.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 8;
                elseif 3 < a_b(i)
                    Ks(i) = 7.5;
                end
    
            elseif b_rt(i) < 3
    
                if  a_b(i) < 1.25
                    Ks(i) = 8.5;
                elseif 1.25 < a_b(i) && a_b(i) < 1.5
                    Ks(i) = 7.5;
                elseif 1.5 < a_b(i) && a_b(i) < 2
                    Ks(i) = 6.5;
                elseif 2 < a_b(i) && a_b(i) < 3
                    Ks(i) = 6;
                elseif 3 < a_b(i)
                    Ks(i) = 5.5;
                end
            else
                disp('b/rt out of bounds')
            end
        
        elseif b_curve(i) > a
            b_a(i) = b_curve(i)./a;
            a_rt(i) = a ./ sqrt(r_DCell(i).*t_DCell);
            
    
            if 8 < a_rt(i) && a_rt(i) < 9
    
                if  b_a(i) < 1.5
                    Ks(i) = 23;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 20;
                else
                    Ks(i) = 18.5;
                end
    
            elseif 7 < a_rt(i) && a_rt(i) < 8
    
                if  b_a(i) < 1.5
                    Ks(i) = 20;
                elseif 1.5 < b_a(i) && b_a(i) < 5
                    Ks(i) = 17.5;
                else
                    Ks(i) = 16;
                end
    
            else
                disp('a/rt out of bounds')
            end
        end
    end
end


% =========================================================================

