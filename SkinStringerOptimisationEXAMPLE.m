clear;
clc;

% Setup material properties and structural constraints
aluminium_density = 2.85e3; % kg/m^3
M_yy_max = 2435955.144; % Maximum bending moment (N*m)
direct_yield_stress = 4.31e8; % Pa
shear_yield_stress = 3.1e8; % Pa
frame_spacing = 0.516; % m
K = 3.62; % Buckling coefficient
aluminium_youngs_modulus = 7.3e10; % Pa

% Optimizer settings (stratified starting values to avoid local minima)
stringer_thickness = linspace(0.00101, 0.0039, 4);
stringer_height = linspace(0.0101, 0.099, 4);
stringer_width = linspace(0.0101, 0.099, 4);
stringer_number = [11, 20, 50, 99]; 
skin_thickness = 0.001;

% Initialize storage for best solutions
x_stress_holder = [];
x_string_holder = [];
x_skin_holder = [];
fval_stress_best = Inf; % Set high initial values
fval_string_best = Inf;
fval_skin_best = Inf;

% Loop through different start points
for i = 1:length(stringer_thickness)
    for j = 1:length(stringer_height)
        for k = 1:length(stringer_width)
            for l = 1:length(stringer_number)
                
                % Initial guess for optimization
                x0 = [stringer_thickness(i), stringer_height(j), stringer_width(k), stringer_number(l)];
                
                % Objective function (minimizing mass)
                mass_func = @(x) ((x(3) * x(1) * 2 + ((x(2) - 2 * x(1)) * x(1))) * aluminium_density ..2.
                                  * x(4) + (2.5512 * pi * skin_thickness * aluminium_density));

                % Stress due to bending constraint
                max_stress_func = @(x) deal([], direct_yield_stress - ...
                    (M_yy_max * (2.5512 / 2) / ...
                    (sum(((2.5512 / 2) * sind(0:360 / x(4):360)).^2 * (2.5512 * pi) / x(4) * skin_thickness + ...
                    (x(3) * x(1) * 2 + (x(2) - 2 * x(1)) * x(1))))));

                % Stringer buckling constraint
                stringer_buckling_func = @(x) deal([], ...
                    -(M_yy_max * (2.5512 / 2) / ...
                    (sum(((2.5512 / 2) * sind(0:360 / x(4):360)).^2 * (2.5512 * pi) / x(4) * skin_thickness + ...
                    (x(3) * x(1) * 2 + (x(2) - 2 * x(1)) * x(1))))) ...
                    - (pi^2 * aluminium_youngs_modulus * ((2 * ((x(2) - x(1)) / 2)^2) * (x(3) * x(1)) + ...
                    x(1) / 12 * (x(2) - 2 * x(1))^3)) / ...
                    ((x(3) * x(1) * 2 + ((x(2) - 2 * x(1)) * x(1))) * frame_spacing^2));

                % Skin buckling constraint
                skin_buckling_func = @(x) deal([], shear_yield_stress - K * aluminium_youngs_modulus * ...
                    ((skin_thickness / ((2.5512 * pi) / x(4)))^2));

                % Define optimization bounds
                Aeq = [];
                Beq = [];
                lb = [0.001, 0.01, 0.01, 10];
                ub = [0.004, 0.1, 0.1, 100];

                % Optimization options
                options = optimoptions('fmincon', 'Display', 'iter');

                % Solve for stress constraints
                [x_stress, fval_stress] = fmincon(mass_func, x0, [], [], Aeq, Beq, lb, ub, max_stress_func, options);

                % Solve for stringer buckling constraints
                [x_string_buck, fval_string_buck] = fmincon(mass_func, x0, [], [], Aeq, Beq, lb, ub, stringer_buckling_func, options);

                % Solve for skin buckling constraints
                [x_skin_buck, fval_skin_buck] = fmincon(mass_func, x0, [], [], Aeq, Beq, lb, ub, skin_buckling_func, options);

                % Check and update best solutions
                if fval_stress < fval_stress_best
                    fval_stress_best = fval_stress;
                    x_stress_holder = x_stress;
                end

                if fval_string_buck < fval_string_best
                    fval_string_best = fval_string_buck;
                    x_string_holder = x_string_buck;
                end

                if fval_skin_buck < fval_skin_best
                    fval_skin_best = fval_skin_buck;
                    x_skin_holder = x_skin_buck;
                end
            end
        end
    end
end

% Display results
disp('Optimal Solution to Resist Stress:');
disp(x_stress_holder);
disp('Optimal Mass Value:');
disp(fval_stress_best);

disp('Optimal Solution to Prevent Stringer Buckling:');
disp(x_string_holder);
disp('Optimal Mass Value:');
disp(fval_string_best);

disp('Optimal Solution to Prevent Skin Buckling:');
disp(x_skin_holder);
disp('Optimal Mass Value:');
disp(fval_skin_best);




