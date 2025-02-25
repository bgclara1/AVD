close all;

%{

 to do:
    
    - fix shear flow
    - fix bending
    - fix stringer opt - done 
    - do heavy frames
    - fuselage fatigue?

%}

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
[InertialLoads2] = getInertialDistro2();

% ======== get air load plot ============================

[xDiscr,aero,aero15] = getAirLoad();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           SF PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  get air SF plot 
aeroSF = getAirSF(xDiscr,aero);

% get fusealge SF plot
[fuselageSF,inertialWithReaction] = getFuselageSF(xDiscr,InertialLoads2);
            % you get a bump at 31 that looks diff to reports but our aero
            % load is more disproportional so when u add it on it makes the
            % thing shoot above zero


[inertialAndAirWithReaction, comboSF] = combineSF(xDiscr,InertialLoads2,aero); % effectively the n=3.75 case


% n = 3.75 SF plot
SF_3_75 = getSF_3_75(inertialAndAirWithReaction) ;   %combines the air and fuselage SF at 3.75 loading

% n = -1.5 SF plot
SF_1_5 = getSF_1_5(xDiscr, aero15,inertialAndAirWithReaction);

% OEI SF plot 
SF_OEI = getOEI(inertialAndAirWithReaction);

% Landing SF plot 
SFland = getLandingSF(xDiscr, inertialAndAirWithReaction);

%Total SF plot
%plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI,SFland);


%------------- Plots --------------------------------
                    % 
                    % 
                    % figure;
                    % plot(xDiscr, InertialLoads, 'LineWidth', 1.5)
                    % title('Inertial Loads')                  
                    % 
                    %  figure;
                    % bar(xDiscr, InertialLoads2)
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
                    % stairs(xDiscr, SF_OEI)
                    % title('OEI SF')
                    % 
                    % 
                    % figure;
                    % stairs(xDiscr, comboSF)
                    % title('combo SF')
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                      STRINGER SKIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------- material/stringer properties ---------------------------

            %stringer
            tensileYieldStressStringer = 490*1e+06;
            shearYieldStressStringer = 310*1e6;
            EStringer = 71*1e9;  

            density = 2780;
            L_eff = 0.5;         
            Ks = 5;       
            nu = 0.33;        
            b=0.02;

            %skin
            tensileYieldStressSkin = 331*1e6;
            shearYieldStressSkin = 290*1e6;
            ESkin = 73e9;      
         

            E = 71e9;
            shearYieldStress = 310*1e6;
            tensileYieldStress = 490*1e+06;

    % REMEMBER TO COPY INTO OPTIMISATION THING

% -----------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                     STRINGER OPTIMISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%[skinThickness, stringerSpacing,h, stringerThickness] = stringerOptimisation(density);
skinThickness = 0.0042;
stringerSpacing = 0.1;
h = 0.05;
stringerThickness = 0.0025;

% Stringer Diagram
[x, y, numStringers] = stringerPlot(stringerSpacing);

%[~, ~, ~, minSkinThickness, ~] = plotShearFlow(shearYieldStressSkin);


% Shear Flow plot
[theta, PlotSFlow, circle, ~, totalSFlow] = plotShearFlow(shearYieldStress);

                figure;
                markerSize = 7;
                scatter(x, y,markerSize, 'filled');
                grid on;
                axis equal; 
                margin = 0.3 * max(abs([x, y])); 
                xlim([-max(abs(x)) - margin, max(abs(x)) + margin]);
                ylim([-max(abs(y)) - margin, max(abs(y)) + margin]);
                xlabel('X-axis', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Y-axis', 'FontSize', 12, 'FontWeight', 'bold');
                title('Stringer Positions Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');


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


stringerArea = ZStringerArea(stringerThickness, h, L); % h height, L flange length

% Direct stress plot
directStress = stressPlot(skinThickness,stringerArea,y,numStringers,density);


                figure;
                hold on;
                grid on;
                plot3(x, zeros(size(x)), y, 'r-', 'LineWidth', 2);
                quiver3(x, zeros(size(x)), y,zeros(size(x)), directStress, zeros(size(x)),'b', 'LineWidth', 1, 'MaxHeadSize', 0.5);
                xlabel('z');
                ylabel('y');
                zlabel('Direct stress \sigma_z (MPa)');
                title('Direct stress distribution around fuselage');
                view(3);
                legend('Fuselage cross-section', 'Direct stress at each stringer');
                axis equal;


%skinThickness = 0.01; %important var for skin bay buckling

%structural compliance
[stringerStressCompliant,stringerEulerCompliant] = checkStringer(directStress, tensileYieldStressStringer,EStringer, b, h, stringerArea, L_eff);
[skinYieldCompliant, skinBucklingCompliant] = checkSkinBay(totalSFlow, skinThickness, stringerSpacing, shearYieldStressSkin, Ks, ESkin, nu);
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

% figure;
% surf(bGrid, hGrid, frameThicknessMatrix', 'EdgeColor', 'none'); % Note the transpose
% xlabel('Flange Width b (m)');
% ylabel('Web Height h (m)');
% zlabel('Frame Thickness t_f (m)');
% title('Frame Thickness with Constant I_{xx}');
% colorbar;
% grid on;
% 
% 
% % 3D Plot for frame area
% figure;
% surf(bGrid, hGrid, frameAreaMatrix', 'EdgeColor', 'none');
% xlabel('Flange Width b (m)');
% ylabel('Web Height h (m)');
% zlabel('Frame Area A (m^2)');
% title('Frame Area with Constant I_{xx}');
% colorbar;
% grid on;

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


angle = pi*(0:5:360)/180;
theta= 1.138;


% Wing front Spar 
angle = linspace(0, 2*pi, 73);
Q = (SF_3_75(31))*cos(theta);  
P = (SF_3_75(31))*sin(theta);           
%T = BM_3_75(31); 
T = 0;

N = Q/(2*pi) * ((pi - angle).*cos(angle) - 0.5*sin(angle)); 
S = (Q/(2*pi)) * (1 + 0.5*cos(angle) - (pi - angle).*sin(angle)) + (P/(2*pi)); 
M = (Q*D)/(4*pi) * (pi - angle) .* (1 - cos(angle)) - 1.5*sin(angle) + (T/(2*pi));

figure;
plot(angle, N);
hold on;
plot(angle, S);
hold on;
plot(angle, M);
xlabel('Angle around frame (rad)');
ylabel('Load Distribution');
legend('Normal Force (N)', 'Shear Force (S)', 'Moment (M)');
title('Front Spar Frame');
grid on;

% Rear front Spar 
angle = linspace(0, 2*pi, 73);
Q = 0;  
P = abs(SF_3_75(39));        
T = BM_3_75(39);        

N = Q/(2*pi) * ((pi - angle).*cos(angle) - 0.5*sin(angle)); 
S = (Q/(2*pi)) * (1 + 0.5*cos(angle) - (pi - angle).*sin(angle)) + (P/(2*pi)); 
M = (Q*D)/(4*pi) * (pi - angle) .* (1 - cos(angle)) - 1.5*sin(angle) + (T/(2*pi));

% figure;
% plot(angle, N);
% hold on;
% plot(angle, S);
% hold on;
% plot(angle, M);
% xlabel('Angle around frame (rad)');
% ylabel('Load Distribution');
% legend('Normal Force (N)', 'Shear Force (S)', 'Moment (M)');
% title('Rear Spar Frame');
% grid on;

phi = pi*(0:5:360)/180; % Angle in radians from 0 to 2*pi
Pt = abs(SF_3_75(31));
Qt = 0; 
Tt = 0;   
R = 3.2425;   

[Np, Sp, Mp, Nq, Sq, Mq, Nt, St, Mt, Nf, Sf, Mf] = calculateFrameLoads(Pt, Qt, Tt, R, phi);

% figure;
% plot(phi, Nf, phi, Sf, phi, Mf);
% legend('Normal Force N_f', 'Shear Force S_f', 'Moment M_f');
% xlabel('Angle \phi (rad)');
% ylabel('Force/Moment');
% title('Combined Frame Loads');
% grid on;

Pt = abs(SF_3_75(39));
Qt = 0; 
Tt = 0;   
R = 3.2425;   

[Np, Sp, Mp, Nq, Sq, Mq, Nt, St, Mt, Nf, Sf, Mf] = calculateFrameLoads(Pt, Qt, Tt, R, phi);

% figure;
% plot(phi, Nf, phi, Sf, phi, Mf);
% legend('Normal Force N_f', 'Shear Force S_f', 'Moment M_f');
% xlabel('Angle \phi (rad)');
% ylabel('Force/Moment');
% title('Combined Frame Loads');
% grid on;
% 







