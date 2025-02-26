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

[xDiscr,aero,aero15,LHT,LHT15,LHTLanding] = getAirLoad();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           SF PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  get air SF plot 
aeroSF = getAirSF(xDiscr,aero);

% get fusealge SF plot
[fuselageSF,inertialWithReaction] = getFuselageSF(xDiscr,InertialLoads);

[inertialAndAirWithReaction375, comboSF375] = combineSF375(xDiscr,InertialLoads,LHT); % effectively the n=3.75 case

[inertialAndAirWithReactionNeg15, comboSFneg15] = combineSFneg15(xDiscr,InertialLoads,aero); % effectively the n=3.75 case

[inertialAndAirWithReactionLanding, comboSFLanding] = combineSFLanding(xDiscr,InertialLoads,LHTLanding); % effectively the n=3.75 case


% n = 3.75 SF plot
SF_3_75 = getSF_3_75(inertialAndAirWithReaction375) ;   %combines the air and fuselage SF at 3.75 loading

% n = -1.5 SF plot
SF_1_5 = getSF_1_5(xDiscr, aero15,inertialAndAirWithReactionNeg15);

% OEI SF plot 
%SF_OEI = getOEI(inertialAndAirWithReaction);

%change param to inertial Load to recalc air for landing


% Landing SF plot 
SFland = getLandingSF(xDiscr, inertialAndAirWithReactionLanding);

%Total SF plot
%plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI,SFland);


%------------- Plots --------------------------------


                    figure;
                    plot(xDiscr, InertialLoads, 'LineWidth', 1.5)
                    title('Inertial Loads')                  

                     figure;
                    bar(xDiscr, InertialLoads2)
                    title('Inertial Loads')

                    figure;
                    bar(xDiscr, aero)
                    hold on;
                    bar(xDiscr,aero15)
                    title('Aero Loads')

                    figure;
                    plot(xDiscr, aeroSF)
                    title('Aero Loads')

                    % figure;
                    % plot(xDiscr, fuselageSF)
                    % title('Fuselage Only SF')

                    figure;
                    plot(xDiscr, SF_3_75)
                    title('SF at load factor 3.75')

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
                    figure;
                    plot(xDiscr, SFland)
                    title('SF at landing')
                    % 
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
    

 %[skinThickness,stringerThickness,h,L,stringerSpacing] = stringerOptimisation(density)


skinThickness = 0.0042;
stringerSpacing = 0.1;
t_s = stringerSpacing;
h = 0.003;
stringerThickness = 0.0025;
L = 0.01;

stringerArea = ZStringerArea(stringerThickness, h, L); % h height, L flange length

% Stringer Diagram
[x, y, numStringers] = stringerPlot(stringerSpacing);

%[~, ~, ~, minSkinThickness, ~] = plotShearFlow(shearYieldStressSkin);


% Shear Flow plot
[theta, PlotSFlow, circle, ~, totalSFlow] = plotShearFlow(shearYieldStress);

                % figure;
                % markerSize = 7;
                % scatter(x, y,markerSize, 'filled');
                % grid on;
                % axis equal; 
                % margin = 0.3 * max(abs([x, y])); 
                % xlim([-max(abs(x)) - margin, max(abs(x)) + margin]);
                % ylim([-max(abs(y)) - margin, max(abs(y)) + margin]);
                % xlabel('X-axis', 'FontSize', 12, 'FontWeight', 'bold');
                % ylabel('Y-axis', 'FontSize', 12, 'FontWeight', 'bold');
                % title('Stringer Positions Around Fuselage', 'FontSize', 14, 'FontWeight', 'bold');
                % 
                % 
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
                % 
                % 


% Direct stress plot
directStress = stressPlot(skinThickness,stringerArea,y,numStringers,density);


                % figure;
                % hold on;
                % grid on;
                % plot3(x, zeros(size(x)), y, 'r-', 'LineWidth', 2);
                % quiver3(x, zeros(size(x)), y,zeros(size(x)), directStress, zeros(size(x)),'b', 'LineWidth', 1, 'MaxHeadSize', 0.5);
                % xlabel('z');
                % ylabel('y');
                % zlabel('Direct stress \sigma_z (MPa)');
                % title('Direct stress distribution around fuselage');
                % view(3);
                % legend('Fuselage cross-section', 'Direct stress at each stringer');
                % axis equal;
                % 

%skinThickness = 0.01; %important var for skin bay buckling



[stringerStressCompliant,stringerEulerCompliant] = checkStringer(directStress, tensileYieldStress,E, h,stringerArea, L, t_s);
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

% [bGrid, hGrid] = meshgrid(bRange, hRange);
% 
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

%   Q = radial
%   P = tangential
%   T = torque

%   N = normal force
%   S = shear force
%   M = moment

% Wing front Spar 

theta = 1.138; %dihedral

angle = pi*(0:5:360)/180;
Q = (SF_3_75(31))*cos(theta);  
P = (SF_3_75(31))*sin(theta);           
T = 0;

[Nf, Sf, Mf] = calculateFrameLoads(P, Q, T, r, angle);
maxNFWS = max(Nf) % maximum N front wing spar
maxSFWS = max(Sf)
maxMFWS = max(Mf)


% 
% figure;
% plot(angle, Nf);
% hold on;
% plot(angle, Sf);
% hold on;
% plot(angle, Mf);
% xlabel('Angle around frame (rad)');
% ylabel('Load Distribution');
% legend('Normal Force (N)', 'Shear Force (S)', 'Moment (M)');
% title('Front Wing Spar Frame');
% grid on;

% Rear wing Spar 
angle = pi*(0:5:360)/180;
Q = (SF_3_75(39))*cos(theta);  
P = (SF_3_75(39))*sin(theta);         
T = 0;

[Nf, Sf, Mf] = calculateFrameLoads(P, Q, T, r, angle);

maxNRWS = max(Nf);% maximum N rear wing spar
maxSRWS = max(Sf);
maxMRWS = max(Mf);

% figure;
% plot(angle, Nf);
% hold on;
% plot(angle, Sf);
% hold on;
% plot(angle, Mf);
% xlabel('Angle around frame (rad)');
% ylabel('Load Distribution');
% legend('Normal Force (N)', 'Shear Force (S)', 'Moment (M)');
% title('Rear Wing Spar Frame');
% grid on;


%front emp frame
r = 2.655;

P = 0;
Q = SF_3_75(71); 
T = SF_OEI(71)/2;   

[Nf, Sf, Mf] = calculateFrameLoads(P, Q, T, r, angle);

maxNFTS = max(Nf);% maximum N front tail spar
maxSFTS = max(Sf);
maxMFTS = max(Mf);

% figure;
% plot(angle, Nf);
% hold on;
% plot(angle, Sf);
% hold on;
% plot(angle, Mf);
% legend('Normal Force N_f', 'Shear Force S_f', 'Moment M_f');
% xlabel('Angle \phi (rad)');
% ylabel('Force/Moment');
% title('Front Tailplane Spar Frame');
% grid on;


%Rear tailplane frame

P = 0;
Q = SF_3_75(75); 
T = SF_OEI(75)/2;   
r = 1.805;

[Nf, Sf, Mf] = calculateFrameLoads(P, Q, T, r, angle);

maxNRTS = max(Nf); % maximum N rear tail spar
maxSRTS = max(Sf);
maxMRTS = max(Mf);
% 
% figure;
% plot(angle, Nf);
% hold on;
% plot(angle, Sf);
% hold on;
% plot(angle, Mf);
% legend('Normal Force N_f', 'Shear Force S_f', 'Moment M_f');
% xlabel('Angle \phi (rad)');
% ylabel('Force/Moment');
% title('Rear Tailplane Spar Frame');
% grid on;
% 

% h     web height
% l     flange length
% t_w   web thickness
% t_f   flange thickness


%[h,l,t] = heavyFrameOptimisation(maxN,maxS,maxM);

disp('Iterating Front Spar Wing ')
[hFWS,lFWS,tFWS] = heavyFrameOptimisation(maxNFWS,maxSFWS,maxMFWS);
disp('Iterating Rear Spar Wing ')
[hRWS,lRWS,tRWS] = heavyFrameOptimisation(maxNRWS,maxSRWS,maxMRWS);
disp('Iterating Front Spar Tail ')
[hFTS,lFTS,tFTS] = heavyFrameOptimisation(maxNFTS,maxSFTS,maxMFTS);
disp('Iterating Rear Spar Tail ')
[hRTS,lRTS,tRTS] = heavyFrameOptimisation(maxNRTS,maxSRTS,maxMRTS);
disp('done ')

% [hFWS,lFWS,tFWS] 
% [hRWS,lRWS,tRWS]
% [hFTS,lFTS,tFTS]
% [hRTS,lRTS,tRTS] 















