close all;

% ======== set default params for plotting ================

close all;

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
[fuselageSF,total] = getFuselageSF(xDiscr,inertialDistro);
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

% Total SF plot
%plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI,SFland);


%------------- Plots --------------------------------
       
                    % 
                    % 
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
                    % 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           BENDING PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BMInert = cumtrapz(fuselageSF);
BM_3_75 = cumtrapz(SF_3_75);
BM_OEI = cumtrapz(SF_OEI);
BM_land = cumtrapz(SFland);
% 
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
%                           S S O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------- material/stringer properties ---------------------------

            shearYieldStress = 187.06*1e+06;
            StringerSpacing = 0.5;
            density = 2.78e+03;
            StringerArea = 0.005;
            tensileYieldStress = 324*1e6;
            E = 73e9;              
            L_eff = 0.5;          
            b = 0.02;             
            h = 0.05;             
            A_str = StringerArea;  

% -----------------------------------------------------------------


% Shear Flow plot
[theta, PlotSFlow, circle, skinThickness] = plotShearFlow(shearYieldStress);

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
                    
% Stringer Diagram
[x,y,numStringers] = stringerPlot(StringerSpacing);

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

% Direct stress plot
directStress = stressPlot(skinThickness,StringerArea,y,numStringers,density);
            
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

% structural compliance
stressCompliant = checkYieldStress(directStress,tensileYieldStress); % 1 true, 0 false
eulerCompliant = checkEulerBuckling(directStress, E, b, h, A_str, L_eff);







