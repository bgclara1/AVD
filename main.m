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
       
%{

                    figure;
                    plot(xDiscr, InertialLoads)
                    title('Inertial Loads')
                    
                    figure;
                    bar(xDiscr, aero)
                    hold on;
                    bar(xDiscr,aero15)
                    title('Aero Loads')
                    
                    figure;
                    plot(xDiscr, aeroSF)
                    title('Aero Loads')
                    
                    figure;
                    plot(xDiscr, fuselageSF)
                    title('Fuselage Only SF')
                    
                    figure;
                    plot(xDiscr, SF_3_75)
                    title('Aero SF at load factor 3.75')
                    
                    figure;
                    plot(xDiscr, SF_1_5)
                    title('Aero SF at load factor - 1.5')
                    
                    figure;
                    plot(xDiscr, SF_OEI)
                    title('OEI SF')
                    
                    figure;
                    plot(xDiscr, SFland)
                    title('SF at landing')
                    
                    figure;
                    plot(xDiscr, SF_3_75)
                    %hold on;
                    %plot(xDiscr,SF_1_5)
                    hold on;
                    plot(xDiscr, SF_OEI)
                    hold on;
                    plot(xDiscr, SFland)
                    title('SF Plot')
                    legend('n=3.75','OEI', 'landing')


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                           BENDING PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















