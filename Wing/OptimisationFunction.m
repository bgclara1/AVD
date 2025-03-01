%% Skin Thickness ==========================================================

E = 72e9;               % Youngs Modulus (Pa)
rho = 2710;              % Density of Al (kg/m^3) 
M = M_pos;              % Bending moment (N/m)
c = linspace(c_root,c_tip,length(y_wing));      % Wing chord (m)
b2 = 0.0758.*c;             % Wing box height (m)
width = 0.6.*c;             % Wing box width  (m)

% Skin thickness with no stringers
t_no_stringers = ((M.*width)./(3.62.*E.*b2)).^(1/3);   % Skin thickness with no stringers (m)
Area_no_stringers = width.*t_no_stringers;          % Panel area with no stringers (m)

% Skin thickness with striners
b1 = 0.4;                     % Stringer pitch   (m)
n = width./b1;              % No. stringers       (/)
d_to_h = 0.3;               % width to height ratio (/)
h = 36e-3;                     % Stringer height         (m) 
d = h*d_to_h;               % Stringer width          (m)
t_s = 1e-3;                    % Stringer thickness      (m)
A_s = t_s*h + 2*t_s*d;      % Z-stringer area formula (m^2)
t_stringers = ((M.*b1.^2)./(3.62.*E.*width.*b2)).^(1/3);       % Thickness with stringers (m)
N = M./b2./width;             % Axial load per unit length (N/m)
sigma = N./(t_stringers);  % Critical global buckling stress (N/m^2)
t_eff = t_stringers + A_s./(1000^2.*b1);  % Effective skin thickness (m)
Area_with_stringers = n.*b1.*t_eff;             % Panel area with stringers (m^2)



% UNCOMMENT THESE PLOTS TO COMPARE AREA OF STRINGER VS NO STRINGER
% Plot to compare affect of stringer vs. no stringers
figure()
plot(y_wing,t_no_stringers)
hold on
plot(y_wing,t_stringers)
hold off

% Plot of different sectional area with and without stringers
figure()
plot(y_wing,Area_no_stringers,DisplayName='No Stringers')
hold on
plot(y_wing,Area_with_stringers,DisplayName='Stringers')
legend()
hold off


% Calculate the necessary ratios to apply the Catchpole correction factor
thickness_ratio = t_s./(1000*t_stringers);
area_ratio = A_s./(1000^2*b1.*t_stringers);

% Use the digitised Catchpole diagram to calculate the correction factor
% for the given vector of ratios
[stress_ratio, new_thickness] = catchpole_interpolation(thickness_ratio, area_ratio);

sigma_crit = stress_ratio.*sigma;        % Local critical buckling stress (N/m^2)

[efficiency, ~] = Farrar_interpolation(thickness_ratio, area_ratio);

rib_spacing_skin = ((efficiency./sigma_crit).^2).*N.*E;      % Optimal rib spacing (m) 


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

% Calculate the thickness of the front and rear spars
t_front = (abs((q_fs.*b2)./(Ks*E))).^(1/3);
t_rear = (abs((q_rs.*b2)./(Ks*E))).^(1/3);

% Calculate the shear stress in the spars
tau_front = q_fs./t_front;
tau_rear = q_rs./t_rear;

figure()
plot(y_wing,t_front)
hold on
plot(y_wing,t_rear)
hold off

figure()
plot(y_wing,tau_front,"DisplayName","Front")
hold on
plot(y_wing,tau_rear,"DisplayName","Rear")
legend()
hold off

tau_skin = q_torque ./ t_stringers;             % Shearstress in skin (Pa)

% Set last value of tau equal to that of the penultimate to aoid dividing by 0
tau_skin(end) = tau_skin(end - 1);

tau_tresca = 125e6;                             % TRESCA shear yield stress (Pa)

combination = (tau_skin./tau_tresca).^2 + sigma./sigma_crit;

figure()
plot(y_wing,combination)
hold off


%% RIB DESIGN ==============================================================

I_skin_panel = (width.*(t_eff).^3)./12 + width.*t_eff.*(b2./2).^2;          % Moment of inertia (m^4)

Crush_load = 0.5.*((M.^2).*rib_spacing_skin.*b2.*(t_eff.*width))./(E.*I_skin_panel.^2);       % Crush load (N)

t_rib = 1e-3;           % Rib thickness (m)

sigma_rib = Crush_load./(t_rib.*width);         % Crush stress (Pa)



%% MASS OF DESIGN / CONSTRAINTS ============================================

obj = t_eff + (b2.*t_rib)./rib_spacing_skin;        % Total equivalent thickness (m)

Limit_sigma_compressive = 400e6;      % Compresssive yield strength (Pa)

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
    rib_positions

end

% =========================================================================