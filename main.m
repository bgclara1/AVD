close all;

% ======== set default params for plotting ================

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',15)
set(groot,'defaulttextfontsize',15)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')

% ======== get inertial distrobution ====================

[InertialLoads] = getInertialDistro();

% ======== get air load plot ============================

[xDiscr,aero,aero15] = getAirLoad();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           SF PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  get air SF plot 
aeroSF = getAirSF(xDiscr,aero);

% get fusealge SF plot
[fuselageSF,total] = getFuselageSF(xDiscr,InertialLoads);
            % you get a bump at 31 that looks diff to reports but our aero
            % load is more disproportional so when u add it on it makes the
            % thing shoot above zero

% n = 3.75 SF plot
SF_3_75 = getSF_3_75(xDiscr,aero,total) ;   %combines the air and fuselage SF at 3.75 loading

% n = -1.5 SF plot
SF_1_5 = getSF_1_5(xDiscr, aero15,total);

% OEI SF plot 
SF_OEI = getOEI(xDiscr,total);

% Landing SF plot 
SFland = getLandingSF(xDiscr, total);

%Total SF plot
plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI,SFland);


%------------- Plots --------------------------------
       


                    % figure;
                    % plot(xDiscr, InertialLoads)
                    % title('Inertial Loads')
                    % 
                    % figure;
                    % bar(xDiscr, aero)
                    % hold on;
                    % bar(xDiscr,aero15)
                    % title('Aero Loads')
                    % 
                    % figure;
                    % plot(xDiscr, aeroSF)
                    % title('Aero Loads')
                    % 
                    % figure;
                    % plot(xDiscr, fuselageSF)
                    % title('Fuselage Only SF')
                    % 
                    % figure;
                    % plot(xDiscr, SF_3_75)
                    % title('Aero SF at load factor 3.75')
                    % 
                    % figure;
                    % plot(xDiscr, SF_1_5)
                    % title('Aero SF at load factor - 1.5')
                    % 
                    % figure;
                    % plot(xDiscr, SF_OEI)
                    % title('OEI SF')
                    % 
                    % figure;
                    % plot(xDiscr, SFland)
                    % title('SF at landing')
                    % % 
                    % figure;
                    % plot(xDiscr, SF_3_75)
                    % %hold on;
                    % %plot(xDiscr,SF_1_5)
                    % hold on;
                    % plot(xDiscr, SF_OEI)
                    % hold on;
                    % plot(xDiscr, SFland)
                    % title('SF Plot')
                    % legend('n=3.75','OEI', 'landing')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           BENDING PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BMInert = cumtrapz(fuselageSF);
BM_3_75 = cumtrapz(SF_3_75);
BM_OEI = cumtrapz(SF_OEI);
BM_land = cumtrapz(SFland);

% figure;
%                     plot(xDiscr, BM_3_75)
%                     %hold on;
%                     %plot(xDiscr,BMInert)
%                     hold on;
%                     plot(xDiscr, BM_OEI)
%                     hold on;
%                     plot(xDiscr, BM_land)
%                     title('BM Plot')
%                     legend('n=3.75','OEI', 'landing')
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                      STRINGER SKIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------- material/stringer properties ---------------------------

            shearYieldStress = 187.06*1e+06;
            StringerSpacing = 0.2;
            density = 2780;
            tensileYieldStress = 324*1e6;
            E = 73e9;              
            L_eff = 0.5;          
            stringerThickness = 0.0016;
            b = 0.02;             
            h = 0.0221; 
            L = 0.0166;
            Ks = 5;           %buckling coeff
            nu = 0.33;        %poisson


% -----------------------------------------------------------------

[theta, PlotSFlow, circle, skinThickness, totalSFlow] = plotShearFlow(shearYieldStress);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                   STRINGER OPTIMISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 3.195; % fuselage radius
C = 2*pi*r;
density = 2.78E+03;

%% Initial Guess for [skinThickness, spacing, t_s, h, L]
x0 = [0.008, 0.6, 0.003, 0.055, 0.03];

%% Realistic Bounds for Variables [lower, upper]
lb = [skinThickness*1.05, 0.4, 0.0013, 0.015, 0.01];
ub = [0.01, 1, 0.005, 0.1, 0.07];

%% Optimization Options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp',...
                       'MaxIterations', 1000, 'MaxFunctionEvaluations', 5000);

%% Run Optimization
[x_opt, mass_opt] = fmincon(@(x)objective(x, density, r), x0, [], [], [], [], lb, ub,...
                            @(x)constraints(x, shearYieldStress, tensileYieldStress, E, L_eff, Ks, nu, r), options);

%% Display Optimized Results
disp('Optimized Variables:');
disp(['Skin Thickness: ', num2str(x_opt(1)), ' m']);
disp(['Stringer Spacing: ', num2str(x_opt(2)), ' m']);
disp(['Stringer Thickness: ', num2str(x_opt(3)), ' m']);
disp(['Stringer Height: ', num2str(x_opt(4)), ' m']);
disp(['Stringer Flange Length: ', num2str(x_opt(5)), ' m']);
disp(['Optimized Fuselage Mass per Unit Length: ', num2str(mass_opt), ' kg/m']);


% ----------------------------------------------------------------------------------


            skinThickness = x_opt(1);
            stringerSpacing = x_opt(2);
            stringerThickness  = x_opt(3);
            h = x_opt(4);
            L = x_opt(5);


stringerArea = ZStringerArea(stringerThickness, h, L); % h height, L flange length

% Shear Flow plot
[theta, PlotSFlow, circle, ~, totalSFlow] = plotShearFlow(shearYieldStress);

                % figure;
                % polarplot(theta,PlotSFlow,'b-','LineWidth', 3)
                % hold on;
                % polarplot(theta,circle,'r-', 'LineWidth', 2)
                % ax = gca; 
                % ax.ThetaZeroLocation = 'top';
                % ax.ThetaDir = 'clockwise';
                % title('Shear Flow Distribution Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');
                % grid on;
                % legend('Shear Flow', 'Fuselage', 'Location', 'best');



% Stringer Diagram
[x,y,numStringers] = stringerPlot(StringerSpacing);
                % 
                % figure;
                % scatter(x, y, 'filled');
                % grid on;
                % axis equal; 
                % margin = 0.3 * max(abs([x, y])); 
                % xlim([-max(abs(x)) - margin, max(abs(x)) + margin]);
                % ylim([-max(abs(y)) - margin, max(abs(y)) + margin]);
                % xlabel('X-axis', 'FontSize', 12, 'FontWeight', 'bold');
                % ylabel('Y-axis', 'FontSize', 12, 'FontWeight', 'bold');
                % title('Stringer Positions Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');

% Direct stress plot
directStress = stressPlot(skinThickness,stringerArea,y,numStringers,density);

            
                % figure;
                % hold on;
                % grid on;
                % plot3(x, zeros(size(x)), y, 'r-', 'LineWidth', 2);
                % quiver3(x, zeros(size(x)), y,zeros(size(x)), directStress, zeros(size(x)),'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
                % xlabel('z');
                % ylabel('y');
                % zlabel('Direct stress \sigma_z (MPa)');
                % title('Direct stress distribution around fuselage');
                % view(3);
                % legend('Fuselage cross-section', 'Direct stress at each stringer');
                % axis equal;

%skinThickness = 0.01; %important var for skin bay buckling

%structural compliance
[stringerStressCompliant,stringerEulerCompliant] = checkStringer(directStress, tensileYieldStress,E, b, h, stringerArea, L_eff);
[skinYieldCompliant, skinBucklingCompliant] = checkSkinBay(totalSFlow, skinThickness, StringerSpacing, shearYieldStress, Ks, E, nu);
structurallyCompliant = all([stringerStressCompliant, stringerEulerCompliant, skinYieldCompliant, skinBucklingCompliant]) % 0 false, 1 true



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                      PRESSURE LOAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------- pressure variables -------------------------------

            D = 6.485;
            r = D/2;
            cabinPressure = 50507.26;
            atmosphericPressure = 0.1013*1e6;
            pressureLoad = (atmosphericPressure-cabinPressure)*1e-5;

% -----------------------------------------------------------------

% fuselage component stresses
[hoopStress, longitudinalStress, sphericalStress] = pressureStresses(D, skinThickness,pressureLoad,atmosphericPressure);

% pressure thickness requirements
[skinThicknessPressure,domeThickness] = pressureThicknesses(D,pressureLoad,atmosphericPressure,tensileYieldStress, nu);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                      LIGHT FRAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------- FRAME GEOMETRIC AND MATERIAL PARAMETERS -----------------

                    E = 7.3*1e9;
                    D = 6.485;
                    L = 0.5; %changes?
                    C_f = 6.25*1e-5;
                    h = 0.05; %changes?
                    M_ult = max(abs(BM_3_75));

% ---------- FRAME STRUCTURAL PROPERTIES -----------------------------

                    EI = C_f*M_ult*D^2/L;
                    I = EI/E;
                    I_xx = I;
                    b = I*12/h^3;
                    A = b*h;

                    bRange = 0.005:0.001:0.06;
                    hRange = 0.03:0.001:0.1;

frameThicknessMatrix = frameThickness(I_xx,bRange,hRange);
frameAreaMatrix = frameArea(frameThicknessMatrix,bRange,hRange);

[bGrid, hGrid] = meshgrid(bRange, hRange);

figure;
surf(bGrid, hGrid, frameThicknessMatrix', 'EdgeColor', 'none'); % Note the transpose
xlabel('Flange Width b (m)');
ylabel('Web Height h (m)');
zlabel('Frame Thickness t_f (m)');
title('Frame Thickness with Constant I_{xx}');
colorbar;
grid on;


% 3D Plot for frame area
figure;
surf(bGrid, hGrid, frameAreaMatrix', 'EdgeColor', 'none');
xlabel('Flange Width b (m)');
ylabel('Web Height h (m)');
zlabel('Frame Area A (m^2)');
title('Frame Area with Constant I_{xx}');
colorbar;
grid on;

% optimal light frame dimensions
[minArea, idx] = min(frameAreaMatrix(:));
[row, col] = ind2sub(size(frameAreaMatrix), idx);
optimal_b = bRange(row);
optimal_h = hRange(col);
optimal_thickness = frameThicknessMatrix(row, col);

% CONSIDER DOING BUCKLING ANALYSIS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                      HEAVY FRAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------- FRAME STRUCTURAL PROPERTIES -----------------------------
           
                    b = 0.02;
                    h = 0.1;
                    A = b*h;
                    y_c = h/2;
                    I_f = b*h^3/12;


angle = pi*(0:10:360)/180;

tanLoad = 1.00E+06; %CHANGE

for i = 1:length(angle)
    N(i) = tanLoad/(2*pi) * ((pi-angle(i))*cos(angle(i))-0.5*sin(angle(i)));
    S(i) = (tanLoad/(2*pi))*(1+0.5*cos(angle(i))-(pi-angle(i))*sin(angle(i)));
    M(i) = (tanLoad*D)/(4*pi)*(pi-angle(i))*(1-cos(angle(i)))-1.5*sin(angle(i));
end


figure;
plot(angle,N);
hold on;
plot(angle,S);
hold on;
plot(angle,M);
xlabel('Angle (rad)')
ylabel('Load (N)')
legend('N','S','M')









