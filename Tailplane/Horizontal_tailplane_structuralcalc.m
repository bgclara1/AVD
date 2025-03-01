% Housekeeping ============================================================
clear
clc
close all

% Formatting
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',18)
set(groot,'defaulttextfontsize',18)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')



%Horizontal Tailplane-Parameters

M_cruise=0.76;
a_cruise=299.517;
V_cruise=M_cruise*a_cruise;
rho_cruise=(8.91/23.77)*1.225;

Va=130.31*sqrt(23.77/8.91);
Vd=V_cruise/0.8

CL_HT_MTOW=-0.17;
CL_HT_MZFW=-0.13;


HT_weight=2998/2 ;
HT_taper=0.3;
HT_span=21.5;
HT_semispan=HT_span/2 ;
HT_AR=4;
HT_S=4*(HT_semispan^2)/(HT_AR) ;
HT_croot=(2*HT_S/(HT_span*(1+HT_taper))) ;
HT_ctip=HT_taper*HT_croot;
n_ult1=3.75;
n_ult2=-1.5 ;
dy=0.01;
E=70*(10^9);
G=30*(10^9);
y=0:dy:HT_semispan;

MTOW=330800;
MZFW=0.85*MTOW;

 Lh_VA_MTOW_n1=calculate(Va,n_ult1,MTOW)
 Lh_VA_MTOW_n2=calculate(Va,n_ult2,MTOW)
% 
% 
 Lh_VD_MTOW_n1=calculate(Vd,n_ult1,MTOW)
 Lh_VD_MTOW_n2=calculate(Vd,n_ult2,MTOW)






%--LOADING DISTRIBUTION------
% 
HT_initial = dist_start(HT_weight*9.81,HT_semispan,HT_taper);
HT_final = HT_taper*HT_initial;
HT_Inertial_dist = trapezoidal(0,-HT_initial,HT_semispan,-HT_final,y);

checkweightHT=((0.5*(HT_final+HT_initial)*HT_semispan)/9.81)*2;


Lh_VA_MTOW_n1_loads=[];
Lh_VA_MTOW_n2_loads=[];


Lh_VD_MTOW_n1_loads=[];
Lh_VD_MTOW_n2_loads=[];

for i=1:length(y)


Lh_VA_MTOW_n1_loads(i)= (2*Lh_VA_MTOW_n1/(pi*HT_semispan))*(1-((y(i)/HT_semispan))^2)^0.5;
Lh_VA_MTOW_n2_loads(i)= (2*Lh_VA_MTOW_n2/(pi*HT_semispan))*(1-((y(i)/HT_semispan))^2)^0.5;


Lh_VD_MTOW_n1_loads(i)= (2*Lh_VD_MTOW_n1/(pi*HT_semispan))*(1-((y(i)/HT_semispan))^2)^0.5;
Lh_VD_MTOW_n2_loads(i)= (2*Lh_VD_MTOW_n2/(pi*HT_semispan))*(1-((y(i)/HT_semispan))^2)^0.5;
end 




figure()
plot(y,HT_Inertial_dist)
hold on
plot(y,Lh_VA_MTOW_n1_loads)
hold on
plot(y,Lh_VA_MTOW_n2_loads)
hold on
plot(y,Lh_VD_MTOW_n1_loads)
hold on
plot(y,Lh_VD_MTOW_n2_loads)
lgd=legend('Inertia','VA n=3.75','VA n=-1.5', 'VD n=3.75','VD n=-1.5');
fontsize(lgd,10,'points')
title('Load Cases')
grid on
grid minor
xlabel('Semi-span (m)','FontSize',13)
ylabel('N/m','FontSize',13)
hold off


%----SHEAR and Bending Distribution--------

totaldist=Lh_VD_MTOW_n2_loads+HT_Inertial_dist;
yreversed=flip(y);
totaldistreversed=flip(totaldist);

shear_force_reversed = cumtrapz(yreversed,totaldistreversed);
bending_moment_reversed= cumtrapz(yreversed,shear_force_reversed);

shear_force = flip(shear_force_reversed);
bending_moment = flip(bending_moment_reversed);

 figure()
 plot(y,-shear_force)
 title('Shear Force')
 grid on
 grid minor
 xlabel('Semi Span')
 ylabel('Shear Force N/m')
 hold off


 figure()
  plot(y,bending_moment)
  title('Bending Moment')
  grid on
  grid minor
   xlabel('Semi Span')
 ylabel('Bending Moment N')
   hold off

   


%----Torque Distributions-----


HT_chord = trapezoidal(0,HT_croot,HT_semispan,HT_ctip,y);

HT_flexuralaxis_location=zeros(1,length(y));

for i=1:length(HT_chord)
    HT_flexuralaxis_location(i)=(0.2 +(0.5*0.5))*HT_chord(i);
end 

cgchords=[];

for i=1:length(HT_chord)
    cgchords(i)=calculate_cg(HT_chord(i));
end 
% 

Torque=zeros(1,length(y));

for i=1:length(y)
    Torque(i)=(Lh_VA_MTOW_n1_loads(i)*(HT_flexuralaxis_location(i)-(0.25*HT_chord(i))) - 1*HT_Inertial_dist(i)*(cgchords(i)-HT_flexuralaxis_location(i)));

end 


    t_rev = flip(Torque);
    torque_rev = cumtrapz(flip(y),t_rev);
    T = flip(torque_rev);

figure()
plot(y,(T))
ylabel('Torque')
xlabel('Semi-Span')
title('Torque Loading')
grid on
grid minor

%arclengthDcell(6.25)
% 
%-----WING BOX GEOMETRY----------

% 
NACA0015=[
0.0000     0.00000
0.0125     0.02367
0.0250     0.03268
0.0500     0.04443
0.0750     0.05250
0.1000     0.05853
0.1500     0.06682
0.2000     0.07172
0.2500     0.07427
0.3000     0.07502
0.4000     0.07254
0.5000     0.06617
0.6000     0.05704
0.7000     0.04580
0.8000     0.03279
0.9000     0.01810
0.9500     0.01008
1.0000     0.00158];

%NACA0015(8,2)---Base half thickness at 20% chord
%NACA0015(14,2)----Base half thickness at 70% chord

wingboxwidth=[];
wingboxheight=[];
frontsparheight_h1=[];
rearsparheight_h2=[];
Dcelllength=[];
Dcellarea=[];


for i=1:length(y)
    wingboxwidth(i)=0.5*HT_chord(i);
    wingboxheight(i)=min( HT_chord(i)*2*NACA0015(8,2)    , HT_chord(i)*2*NACA0015(14,2)        );

    frontsparheight_h1(i)=HT_chord(i)*2*NACA0015(8,2);
    rearsparheight_h2(i)=HT_chord(i)*2*NACA0015(14,2);
    Dcelllength(i)=arclengthDcell(HT_chord(i))*2;
    Dcellarea(i)=HT_chord(i)*abs(trapz(NACA0015(1:8,1),NACA0015(1:8,2)))*2;
end 

frontsparx=[0.2,0.2,0.2];
frontspary=[ 0.07172,0,- 0.07172];

rearsparx=[0.7,0.7,0.7];
rearspary=[ 0.04580,0,- 0.04580];

idealisedrectanglex=[0.2,0.7];
idealisedrectangley=[0.07172,0.04580];

figure()
plot(NACA0015(:,1),NACA0015(:,2),'b')
hold on
plot(NACA0015(:,1),-NACA0015(:,2),'b')
hold on
plot(frontsparx,frontspary,'r')
hold on
plot(rearsparx,rearspary,'r')
hold on
plot(idealisedrectanglex,idealisedrectangley,'r','LineStyle','--')
hold on
plot(idealisedrectanglex,-idealisedrectangley,'r','LineStyle','--')
ylim([-0.5 0.5])
hold off

% 
% %-------Skin and Stringer-----------

% bstart=linspace(0.01,0.5,1000);
% 
% 
% stressatroot=[];
% 
% 
% for i=1:length(bstart)
%     tstartatroot=((abs(bending_moment(i))*bstart(i)^2)/(wingboxheight(i)*wingboxwidth(i)*3.62*E))^(1/3);
%     stressatroot(i)=abs(bending_moment(i))/(wingboxwidth(i)*wingboxheight(i)*tstartatroot);
% end 
% 
% figure()
% plot(bstart,stressatroot/(10^6))
% yline(310,'LineWidth',2,'LineStyle','--')
% hold on
% grid on
% grid minor
% title('Stringer Pitch vs Stress MPA','FontSize',14)
% legend('Stress at Root','Material Compressive Strength')
% xlabel('Stringer Pitch-m','FontSize',14)
% ylabel('Stress-MPa','FontSize',14)
% hold off


% 
%AS_BT and ts_t arrays correspond to farrar efficiency of 0.7. Stress ratio
%arrays correspond to above AS_BT and ts_t arrays. 24 values in total to
%iterate over

%original
 AS_btratioarray=[0.858044164	0.788643533	0.719242902	0.65615142	0.599369085	0.536277603	0.498422713	0.429022082	0.378548896	0.479495268	0.719242902	0.990536278	1.20189];
 ts_tratioarray=[1.341992883	1.255516014	1.156227758	1.03772242	0.954448399	0.86797153	0.755871886	0.640569395	0.528469751	0.416370107	0.416370107	0.438790036	0.4788];
 stressratioarray=[1.538209983	1.485542169	1.364716007	1.358519793	1.327538726	1.259380379	1.200516351	1.138554217	1.023924269	0.360929432	0.159552496	0.103786575	0.153356282];


%  AS_btratioarray=[0.741433022	0.710280374	0.697819315	0.679127726	0.654205607	0.629283489	0.598130841	0.579439252	0.554517134	0.542056075	0.510903427	0.492211838	0.46728972	0.448598131	0.429906542	0.411214953	0.386292835	0.367601246	0.342679128	0.317757009	0.299065421	0.280373832	0.267912773	0.267912773	0.858044164	0.788643533	0.719242902	0.65615142	0.599369085	0.536277603	0.498422713	0.429022082	0.378548896	1.07082153	1.031161473	0.985835694	0.92917847	0.889518414	0.83286119	0.781869688	0.742209632	0.696883853	0.657223796	0.617563739	0.57223796	0.532577904	0.504249292	0.464589235
% ];
%   ts_tratioarray=[1.411327434	1.369911504	1.325309735	1.280707965	1.245663717	1.210619469	1.162831858	1.121415929	1.070442478	1.03539823	0.993982301	0.94619469	0.90159292	0.850619469	0.799646018	0.761415929	0.710442478	0.665840708	0.624424779	0.58300885	0.528849558	0.487433628	0.446017699	0.40460177	1.341992883	1.255516014	1.156227758	1.03772242	0.954448399	0.86797153	0.755871886	0.640569395	0.528469751	1.38	1.326315789	1.288421053	1.237894737	1.2	1.130526316	1.064210526	1.023157895	0.966315789	0.896842105	0.833684211	0.770526316	0.713684211	0.66	0.581052632
% ];
%   stressratioarray=[1.569473684	1.541052632	1.515789474	1.503157895	1.496842105	1.490526316	1.449473684	1.427368421	1.398947368	1.389473684	1.370526316	1.338947368	1.304210526	1.272631579	1.244210526	1.218947368	1.187368421	1.152631579	1.127368421	1.098947368	1.064210526	1.035789474	1.013684211	0.855789474	1.538209983	1.485542169	1.364716007	1.358519793	1.327538726	1.259380379	1.200516351	1.138554217	1.023924269	1.528209765	1.531464738	1.482640145	1.505424955	1.463110307	1.427305606	1.358951175	1.355696203	1.332911392	1.274321881	1.179927667	1.16039783	1.150632911	1.085533454	1.072513562
% ];



biterate=linspace(0.145,0.5,length(ts_tratioarray));
matrixfordifferentcombo=zeros(length(biterate),length(AS_btratioarray));
combinationumber=linspace(1,length(stressratioarray),length(stressratioarray));

for j=1:length(biterate)
for i=1:length(stressratioarray)
    [panelareanostringers,panelareastringers_1,L_1,effectivethicknesst,thicknessnostringers,ribV1,skinV1,t,stress,AS]=stringerskin(i,biterate(j),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
    matrixfordifferentcombo(j,i)=(ribV1+ skinV1)*2780;
end 
end 

figure()
surf(combinationumber,biterate,matrixfordifferentcombo)
xlabel('Combination Number','FontSize',12);
ylabel('Stringer Pitch-m','FontSize',12);
zlabel('Skin and Rib Weight-kg','FontSize',12);
grid on
grid minor
colormap('jet');
colorbar;
title('Optimisation 1','FontSize',15)
hold off
% 
% 
% 
%trial
 AS_btratioarray=[0.858044164	0.788643533	0.719242902	0.65615142	0.599369085	0.536277603	0.498422713	0.429022082	0.378548896];
 ts_tratioarray=[1.341992883	1.255516014	1.156227758	1.03772242	0.954448399	0.86797153	0.755871886	0.640569395	0.528469751];
 stressratioarray=[1.538209983	1.485542169	1.364716007	1.358519793	1.327538726	1.259380379	1.200516351	1.138554217 1.023924269];

%[panelareanostringers1,panelareastringers1,L1,effectivethicknesst,criticalshearbuckling,thicknessnostringers,ribV,skinV]=stringerskin(index,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);

% figure()
% plot(y,L1)
% ylim([0 5])
% 
%  figure()
%  plot(y,panelareastringers1)
%  hold on
%  plot(y,panelareanostringers1)
%  title('Panel Area with and without Stringers')
%  legend('With Stringers','No Stringers')
%  hold off



 % PANELAREANOSTRINGERS=zeros(length(y),length(stressratioarray));
 % PANELAREASTRINGERS=zeros(length(y),length(stressratioarray));
 % LRESULT=zeros(length(y),length(stressratioarray));

 %iterate to find best combo for a given b


biterate=linspace(0.145,0.5,length(ts_tratioarray));
matrixfordifferentcombo=zeros(length(biterate),length(AS_btratioarray));
combinationumber=[1 2 3 4 5 6 7 8 9];

for j=1:length(biterate)
for i=1:length(stressratioarray)
    [panelareanostringers,panelareastringers_1,L_1,effectivethicknesst,thicknessnostringers,ribV1,skinV1,t,stress,AS]=stringerskin(i,biterate(j),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
    matrixfordifferentcombo(j,i)=(skinV1+ribV1)*2780;
end 
end 

figure()
surf(combinationumber,biterate,matrixfordifferentcombo)
xlabel('Combination Number','FontSize',12);
ylabel('Stringer Pitch-m','FontSize',12);
zlabel('Skin and Rib Weight-kg','FontSize',12);
grid on
grid minor
colormap('jet');
colorbar;
title('Optimisation 2','FontSize',15)
hold off

% b=0.1;
% 
% [panelareanostringers,panelareastringers_1,L_1,effectivethicknesst,thicknessnostringers,ribV1,skinV1,t,stress,AS]=stringerskin(1,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_2,L_2,effectivethicknesst,thicknessnostringers,ribV2,skinV2,t,stress,AS]=stringerskin(2,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_3,L_3,effectivethicknesst,thicknessnostringers,ribV3,skinV3,t,stress,AS]=stringerskin(3,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_4,L_4,effectivethicknesst,thicknessnostringers,ribV4,skinV4,t,stress,AS]=stringerskin(4,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_5,L_5,effectivethicknesst,thicknessnostringers,ribV5,skinV5,t,stress,AS]=stringerskin(5,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_6,L_6,effectivethicknesst,thicknessnostringers,ribV6,skinV6,t,stress,AS]=stringerskin(6,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_7,L_7,effectivethicknesst,thicknessnostringers,ribV7,skinV7,t,stress,AS]=stringerskin(7,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_8,L_8,effectivethicknesst,thicknessnostringers,ribV8,skinV8,t,stress,AS]=stringerskin(8,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% [panelareanostringers,panelareastringers_9,L_9,effectivethicknesst,thicknessnostringers,ribV9,skinV9,t,stress,AS]=stringerskin(9,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % [panelareanostringers,panelareastringers_10,L_10,effectivethicknesst,thicknessnostringers,ribV10,skinV10,t,stress,AS]=stringerskin(10,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % [panelareanostringers,panelareastringers_11,L_11,effectivethicknesst,thicknessnostringers,ribV11,skinV11,t,stress,AS]=stringerskin(11,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % [panelareanostringers,panelareastringers_12,L_12,effectivethicknesst,thicknessnostringers,ribV12,skinV12,t,stress,AS]=stringerskin(12,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % [panelareanostringers,panelareastringers_13,L_13,effectivethicknesst,thicknessnostringers,ribV13,skinV13,t,stress,AS]=stringerskin(13,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% 
% 
% 
% 
% % combination=[1,2,3,4,5,6,7,8,9,10,11,12,13];
% % ribvresults=[ribV1,ribV2,ribV3,ribV4,ribV5,ribV6,ribV7,ribV8,ribV9,ribV10,ribV11,ribV12,ribV13];
% % skivresults=[skinV1,skinV2,skinV3,skinV4,skinV5,skinV6,skinV7,skinV8,skinV9,skinV10,skinV11,skinV12,skinV13];
% % 
% % figure()
% % plot(combination,ribvresults)
% % hold on
% % plot(combination,skivresults)
% % hold on
% % plot(combination,ribvresults+skivresults)
% % xlim([1 13])
% % legend('Rib Volume','Skin Stringer Volume','Total Volume')
% % xlabel('Combination','FontSize',14)
% % ylabel('Volume m3','FontSize',14)
% % title('Optimum Combination of AS/bt,ts/t and Stress Ratio','FontSize',11)
% % grid on
% % grid minor
% % hold off
% % 
% % % 
% % % %iterate to find best b for a given combo
% % % 
% %  biterate=[0.1 0.2 0.3 0.4 0.5 0.6 0.7];
% % % 
% %  [panelareanostringers,panelareastringers_24_1,L_24_1,effectivethicknesst,thicknessnostringers,ribV1_2,skinV1_2,t,stress]=stringerskin(8,biterate(1),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% %  [panelareanostringers,panelareastringers_24_2,L_24_2,effectivethicknesst,thicknessnostringers,ribV2,skinV2,t,stress]=stringerskin(8,biterate(2),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% % % 
% %  [panelareanostringers,panelareastringers_24_3,L_24_3,effectivethicknesst,thicknessnostringers,ribV3,skinV3,t,stress,AS]=stringerskin(8,biterate(3),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% % % 
% %  [panelareanostringers,panelareastringers_24_4,L_24_4,effectivethicknesst,thicknessnostringers,ribV4,skinV4,t,stress,AS]=stringerskin(8,biterate(4),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% % % 
% %  [panelareanostringers,panelareastringers_24_5,L_24_5,effectivethicknesst,thicknessnostringers,ribV5,skinV5,t,stress,AS]=stringerskin(8,biterate(5),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% %  [panelareanostringers,panelareastringers_24_6,L_24_6,effectivethicknesst,thicknessnostringers,ribV6,skinV6,t,stress,AS]=stringerskin(8,biterate(6),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % % 
% %  [panelareanostringers,panelareastringers_24_7,L_24_7,effectivethicknesst,thicknessnostringers,ribV7,skinV7,t,stress,AS]=stringerskin(8,biterate(7),AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);
% % 
% %  ribvresults=[ribV1_2,ribV2,ribV3,ribV4,ribV5,ribV6,ribV7];
% %  skinstringerresults=[skinV1_2,skinV2,skinV3,skinV4,skinV5,skinV6,skinV7];
% % 
% % figure()
% % plot(biterate,ribvresults)
% % hold on
% % plot(biterate,skinstringerresults)
% % hold on
% % plot(biterate,skinstringerresults+ribvresults)
% % xlim([ biterate(1) 0.7])
% % xlabel('Stringer Pitch-m','FontSize',14)
% % ylabel('Volume m3','FontSize',14)
% % title('Optimum Value of Stringer Pitch','FontSize',11)
% % legend('Rib Volume','Skin Stringer Volume','Total Volume')
% % grid on
% % grid minor
% % hold off
% % 
% % 
% % % 
  [panelareanostringers,panelareastringers_24_7,Lfinal,effectivethicknesst,thicknessnostringers,ribV7,skinV7,t,stress,AS]=stringerskin(8,0.2,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth);

ts=[];

for i=1:length(y)
    ts(i)=ts_tratioarray(8)*t(i) ;
end 

tsactual=t*ts_tratioarray(8)*1000 ;
tsactual(length(tsactual))=1;
tsmanufacturable=ceil(t*ts_tratioarray(8)*1000 );
tsmanufacturable(length(tsmanufacturable))=1;

figure()
 plot(y,t*1000)
 hold on
 plot(y,ceil(t*1000))
 hold on
plot(y,tsactual)
hold on
plot(y,tsmanufacturable)
grid on
grid minor
ylabel('Thickness-mm','FontSize',14)
xlabel('Semi-Span m','FontSize',14)
lgd=legend('T theoretical','T Actual','Ts theoretical','Ts Actual');
fontsize(lgd,15,'points')
title('Thickness Distributions-mm','FontSize',14)
hold off
% 

 figure()
 plot(y,Lfinal)
 xlabel('Semi-Span m','FontSize',17)
 grid on
 grid minor
 ylabel('Spacing - m','FontSize',17)
 title('Rib Spacing-m')
 ylim([ 0.6 1])


steppedthickness=generate_step_inputs(y,t);

figure()
plot(y,t*1000)
hold on
plot(y,ceil(t*1000))
hold off



z=ceil(t*1000);





z(length(z))=1;
tfinal=z/1000;

figure()
plot(y,t*1000)
hold on
plot(y,tfinal*1000)
grid on
grid minor
xlabel('Semi-Span (m)','FontSize',17)
ylabel('Thickness mm','FontSize',17)
title('Skin Thickness Distribution','FontSize',17)
legend('Theoretical','Actual')
hold off


figure()
plot(y,AS*10000)
hold on
plot(y,ceil(AS*10000))
lgd=legend('Calculated Stringer Area','Actual Stringer Area');
fontsize(lgd,15,'points')
xlabel('Semi-Span m','FontSize',17)
ylabel('Area of Stringers-cm 2','FontSize',17)
title('Stringer Area Distribution','FontSize',17)
grid on 
grid minor
hold off

tfrontcalc=[];
trearcalc=[];
q1calc=[];
q2calc=[];


shearstressfrontspar_actual=[];
shearstressrearspar_actual=[];
% 
% 
for i=1:length(y)-1

torque=abs(T(i));
Fy=abs(shear_force(i)) ;
frontsparheight=frontsparheight_h1(i);
rearsparheight=rearsparheight_h2(i);
wingboxwidthcalc=wingboxwidth(i);
G=G;
thicknessskin=t(i);
Dcell_length=Dcelllength(i);
Dcellarea_calc=Dcellarea(i);



initial_guess = [5e-3; 5e-3;Fy/2*frontsparheight ; Fy/2*rearsparheight];

% Solve the system
options = optimoptions('fsolve', 'Display', 'iter'); % Optional: Display iterations
solution = fsolve(@(variables) Torque2(variables,torque,Fy,frontsparheight,rearsparheight,wingboxwidthcalc,G,E,thicknessskin,Dcell_length,Dcellarea_calc) , initial_guess, options);

tfrontcalc(i)=solution(1);
trearcalc(i)=solution(2);
q1calc(i)=solution(3);
q2calc(i)=solution(4);


shearstressfrontspar_actual(i)=abs(q2calc(i)-q1calc(i)-(Fy/2*frontsparheight))/tfrontcalc(i);
shearstressrearspar_actual(i)=abs(q1calc(i)-(Fy/2*rearsparheight))/trearcalc(i) ;



end 


stress1=[];

for i=1:length(y)
    stress1(i)=(abs(bending_moment(i))/(wingboxwidth(i)* HT_chord(i)*2*((NACA0015(8,2)+ NACA0015(8,2)  )*0.5)))/tfinal(i) ;
end 

figure()
plot(y,stress1/(10^6))
yline(290,'LineWidth',2,'LineStyle','--')
xlabel('Semi Span - m','FontSize',17)
grid on 
grid minor
legend('Applied Stress','Compressive Yield Stress')
title('Compressive Stress Distribution','FontSize',17)
ylabel(' Stress MPa','FontSize',17)
ylim([0 320])


shearstressfrontspar_actual=remove_negatives(shearstressfrontspar_actual);
shearstressrearspar_actual=remove_negatives(shearstressrearspar_actual);


figure()
plot(linspace(0,HT_semispan,length(shearstressfrontspar_actual)),shearstressfrontspar_actual/(10^6))
hold on
plot(linspace(0,HT_semispan,length(shearstressrearspar_actual)),shearstressrearspar_actual/(10^6))
yline(310,'LineWidth',2,'LineStyle','--')
legend('Front Spar','Rear Spar','Yield Stress')
grid on
grid minor
xlabel('Semi-Span (m)','FontSize',17)
ylabel('Stress-MPA','FontSize',17)
title('Shear Stress on Spars','FontSize',17)
ylim([0 320])
hold off
% 
% % %-------RIBS-------------
% 
ribthickness=[];

for i=1:length(y)



    I=((HT_chord(i)*(thicknessnostringers(i)^3))/12) + (HT_chord(i)*t(i)*(wingboxheight(i)/2)^2) ;
    F=(bending_moment(i)^2)*Lfinal(i)*wingboxheight(i)*HT_chord(i)*t(i)/(2*E*I^2) ;
    ribthickness(i)=(F*wingboxheight(i)/(3.62*E))^(1/3);

end 

x = 0; % Starting point
endpoint = HT_semispan; % Ending point
z = Lfinal; % Array of spacing values

% Compute the cumulative sum of the spacing array
cumulativeSpacing = cumsum(z);

% Generate the new array by adding the starting point
ribspacingy = x + [0, cumulativeSpacing];

% Ensure the last point does not exceed the ending point y
ribspacingy (ribspacingy  > endpoint) = []; % Remove points beyond y


% Define x and y arrays
x = y; % x-axis values
yarray = ceil(ribthickness*1000); % y-axis values

% Define specific x points (not necessarily integers)
xPoints = ribspacingy;

% Interpolate y at the specified x points
yPoints = interp1(x, yarray, xPoints, 'linear'); % Linear interpolation

% Plot the points
figure;
plot(xPoints, yPoints, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
title('Plot of y at Specific x Points (Interpolated)');
xlabel('x');
ylabel('y');
grid on;


figure()
plot(y,ribthickness*1000)
hold on
plot(y,ceil(ribthickness*1000))
hold on
plot(xPoints,yPoints,'Marker','x','LineStyle','none','MarkerSize',10)
lgd=legend('Theoretical Rib Thickness','Actual Rib Thickness','Rib Locations');
fontsize(lgd,14,'points')
grid on
grid minor
xlabel('Semi-Span m','FontSize',17)
ylabel('Thickness - mm','FontSize',17)
title('Rib Thickness-mm','FontSize',17)
hold off




stringersperspan=[];

for i=1:length(y)
    stringersperspan(i)=wingboxwidth(i)/0.2 ;
end 

% Define x and y arrays
x = y; % x-axis values
yarray = ceil(stringersperspan); % y-axis values

% Define specific x points (not necessarily integers)
xPoints = ribspacingy;

% Interpolate y at the specified x points
yPoints = interp1(x, yarray, xPoints, 'linear'); % Linear interpolation



figure()
plot(y,ceil(stringersperspan))
hold on
plot(xPoints,yPoints,'Marker','x','LineStyle','none','MarkerSize',10)
grid on
grid minor
legend('Stringers','Rib Locations')
xlabel('Semi Span m')
ylabel('Number of Stringers')
title('Stringers Distribution')
hold off

% 

% %------Spars---------------




criticalshearbucklingactual=[];
ratio1=[];

for i=1:length(y)
    skinshearstress_skin(i)=(abs(T(i))/(0.5*(frontsparheight_h1(i)+rearsparheight_h2(i))*wingboxwidth(i)))/(tfinal(i));
    criticalshearbucklingactual(i)=8.1*E*(tfinal(i)/0.2)^2 ;

   if criticalshearbucklingactual(i)>(300e6)
       criticalstress=300e6 ;
   end 

   if criticalshearbucklingactual(i)<(300e6)
       criticalstress=criticalshearbucklingactual(i);
   end 


    ratio1(i)=((skinshearstress_skin(i)/criticalstress)^2)  + (1/stressratioarray(8)) ;
end 

ratio2=[];

stressratioarraynew=[];
ASarraynew=[];
tsarraynew=[];

for i=1:length(y)
    if y(i) >= 1
        stressratioarraynew(i)=stressratioarray(8) ;
        ASarraynew(i)=AS_btratioarray(8)*0.2*t(i);
        tsarraynew(i)=ts_tratioarray(8)*t(i);
    end 
    if y(i) < 1
        stressratioarraynew(i)=stressratioarray(7) ;
        ASarraynew(i)=AS_btratioarray(7)*0.2*t(i);
        tsarraynew(i)=ts_tratioarray(7)*t(i);
    end 

end 

for i=1:length(y)-1
    skinshearstress_skin_2(i)=q1calc(i)/tfinal(i) ;

   if criticalshearbucklingactual(i)>(300e6)
       criticalstress=300e6 ;
   end 

   if criticalshearbucklingactual(i)<(300e6)
       criticalstress=criticalshearbucklingactual(i);
   end 

    ratio2(i)=((skinshearstress_skin_2(i)/criticalstress)^2)  + (1/stressratioarraynew(i)) ;
end 

RcRs=[];

for i=1:length(y)
    RcRs(i)=0.99 ;
end 

figure()
plot(y,tsarraynew*1000)
hold on
plot(y,ceil(tsarraynew*1000))
grid on
grid minor
title('Stringer Skin Thickness Distribution New')
xlabel('Semi-Span m')
ylabel('Stringer Skin Thickness - mm')
hold off


figure()
plot(y,ASarraynew*10000)
hold on
plot(y,ceil(ASarraynew*10000))
grid on
grid minor
title('Stringer Area Distribution New')
xlabel('Semi-Span m')
ylabel('Stringer Area - cm 2')
hold off



figure()
plot(y,ratio1)
hold on
plot(y(1:end-1),ratio2)
hold on
plot(y,RcRs,'LineStyle','--')
hold off
legend('Without D Cell','With D Cell','Combined Loading Condition')
title('Combined Shear and Compression Loading Condition-New','FontSize',17)
xlabel('Semi-Span (m)','FontSize',17)
ylabel('Combined Loading Ratio','FontSize',17)
grid on
grid minor
hold off

figure()
plot(y(1:end-1),trearcalc*1000)
hold on
plot(y(1:end-1),tfrontcalc*1000)
hold off



trearcalcfiltered=remove_negatives(trearcalc);

trearcalcfilteredd=ceil(trearcalcfiltered*1000);
trearcalcfilteredd(trearcalcfilteredd(1:400) == 4) = 5;

tfrontcalcfiltered=remove_negatives(tfrontcalc);

frontsparheightfilt=linspace(frontsparheight_h1(1),frontsparheight_h1(length(frontsparheight_h1)),length(tfrontcalcfiltered))
rearsparheightfilt=linspace(rearsparheight_h2(1),rearsparheight_h2(length(rearsparheight_h2)),length(trearcalcfilteredd))






 yfiltered=linspace(0,HT_semispan,length(trearcalcfiltered));

areafront=[];
arearear=[];
for i=1:length(yfiltered)
    areafront(i)=frontsparheightfilt(i)*tfrontcalcfiltered(i);
    arearear(i)=rearsparheightfilt(i)*trearcalcfiltered(i);
end 


frontsparweight=abs(trapz(yfiltered,areafront))*2780
rearsparweight=abs(trapz(yfiltered,arearear))*2780



figure()
plot(linspace(0,HT_semispan,length(trearcalcfiltered)),trearcalcfiltered*1000)
hold on
plot(yfiltered,trearcalcfilteredd)
xlabel('Semi-Span m','FontSize',17)
ylabel('Thickness mm','FontSize',17)
title('Rear Spar Thickness Distribution','FontSize',17)
grid on 
grid minor
legend('Calculated','Actual')
hold off


figure()
plot(linspace(0,HT_semispan,length(tfrontcalcfiltered)),tfrontcalcfiltered*1000)
hold on
plot(linspace(0,HT_semispan,length(tfrontcalcfiltered)),ceil(tfrontcalcfiltered*1000))
xlabel('Semi-Span m','FontSize',17)
ylabel('Thickness mm','FontSize',17)
title('Front Spar Thickness Distribution','FontSize',17)
grid on 
grid minor
legend('Calculated','Actual')
hold off



% 
figure()
plot(y(1:end-1),skinshearstress_skin_2/(10^6))
hold on
plot(y,criticalshearbucklingactual/(10^6))
hold on
yline(300,'LineWidth',2,'LineStyle','--')
grid on
grid minor
lgd=legend('Skin Shear Stress','Critical Shear Buckling Stress','Shear Yield Stress');
fontsize(lgd,14,'points')
title('Shear Stresses in Skin')
xlabel('Semi-Span (m)')
ylabel('Stress MPa')
hold off
% 
figure()
plot(y,skinshearstress_skin/(10^6))
hold on
plot(y,criticalshearbucklingactual/(10^6))
legend('Skin Shear Stress','Critical Shear Buckling Stress')
title('Shear Stresses in Skin- MPA')
hold off

% % 
% % 
% %-----D Cell--------------
% 

Dcellstressdistribution=[];

for i=1:length(y)-1
    Dcellstressdistribution(i)=abs(q2calc(i))/(10e-3);
end 

plot(y(1:end-1),Dcellstressdistribution/(10^6))

q2calcfiltered=remove_negatives(q2calc);
ycalcq2=linspace(0,HT_semispan,length(q2calcfiltered));

q2calcdcell=[5.268e05,4.1928e05,3.1829e05,2.3153e05,1.6538e05,1.0787e05,6.1116e04,2.9196e04,8.5518e03,53.6004];

tresca=150e6;
tdcell=5e-3;
psuedoribnumbers=[1,3,5,7,9,11,13,15];

dcellb=[];
b_a=[];

yDcell=linspace(0,HT_semispan,10);

HTchorddcell=trapezoidal(0,HT_croot,HT_semispan,HT_ctip,yDcell);

a_bpsuedorib=zeros(length(HTchorddcell),length(psuedoribnumbers)) ;
b_rt=zeros(length(HTchorddcell),length(psuedoribnumbers)) ;

b_apseudorib=zeros(length(HTchorddcell),length(psuedoribnumbers)) ;
a_rt=zeros(length(HTchorddcell),length(psuedoribnumbers)) ;

%find Ks distirbution--->find thickness distribution----> tradeoff between
%ribs and panel thickness

for j=1:length(psuedoribnumbers)
    a=HT_semispan/(psuedoribnumbers(j)+1);
for i=1:length(yDcell)
b=arclengthDcell(HTchorddcell(i));

R=HTchorddcell(i)*NACA0015(8,1);

if (a/b) > 1
    a_bpsuedorib(i,j)=a/b ;
    b_rt(i,j)=b/sqrt(R*tdcell);
end 

if (a/b) < 1
    b_apseudorib(i,j)=b/a ;
    a_rt(i,j)=a/sqrt(R*tdcell);
end 


end 
end 

KS=[30.60773481	39.44751381	44.69043152	30.84427767	21.05065666	15.87242026	11.81988743	10.01876173
29.28176796	35.13812155	43.31491713	31.2945591	22.28893058	16.88555347	12.49530957	11.14446529
27.9558011	32.26519337	40.99447514	33.0956848	24.76547842	17.89868668	13.50844278	12.1575985
26.29834254	28.83977901	35.5801105	36.13508443	26.34146341	19.9249531	15.30956848	13.84615385
24.19889503	26.85082873	31.04972376	35.13812155	28.93058161	22.06378987	17.44840525	15.19699812
22.76243094	24.30939227	27.0718232	30.82872928	31.2945591	24.09005629	19.69981238	16.43527205
19.33701657	21.32596685	23.42541436	27.9558011	29.61325967	26.79174484	22.06378987	18.12382739
17.56906077	19.1160221	20.77348066	22.43093923	25.3038674	27.40331492	25.10318949	20.71294559
15.80110497	15.80110497	18.12154696	18.78453039	20.11049724	22.09944751	24.53038674	24.54033771
13.59116022	13.14917127	14.58563536	14.91712707	16.57458564	17.45856354	19.55801105	20.11049724
];

dcellpanelvolumes=[];
dcellpanelthicknesses=[];

for j=1:length(psuedoribnumbers)
for i=1:length(yDcell)
        dcellpanelthicknesses(i,j)=((arclengthDcell(HTchorddcell(i))^2)^(1/3))*(q2calcdcell(i)/((KS(i,j)*E)))^(1/3);
        %dcellpanelthicknesses(i,j)=((arclengthDcell(HTchorddcell(i))^2)^(1/2))*(tresca/((KS(i,j)*E)))^(1/2);

        dcellpanelvolumes(i,j)=(arclengthDcell(HTchorddcell(i)))*dcellpanelthicknesses(i,j);
end 
end 

dcellpanelvolumes_differentribnumbers=[];

for i=1:length(psuedoribnumbers)
    dcellpanelvolumes_differentribnumbers(i)=abs(trapz(dcellpanelvolumes(1:end,i),yDcell));
end 


dcellribvolumes=[];

for j=1:length(psuedoribnumbers)
    ydcellribvolumes=linspace(0,HT_croot,psuedoribnumbers(j))
    HTchordRib=trapezoidal(0,HT_croot,HT_semispan,HT_ctip,ydcellribvolumes);
    for i=1:length(ydcellribvolumes)
        dcellribvolumes(i,j)=arclengthDcell(HTchordRib(i))*tdcell ;
    end 

end 

dcellvolumes_differentribnumbers=[];

for i=1:length(psuedoribnumbers)
dcellvolumes_differentribnumbers(i)=sum(dcellribvolumes(1:end,i));
end 

figure()
plot(psuedoribnumbers,dcellvolumes_differentribnumbers*2780)
hold on
plot(psuedoribnumbers,dcellpanelvolumes_differentribnumbers*2780)
hold on
plot(psuedoribnumbers,(dcellpanelvolumes_differentribnumbers+dcellvolumes_differentribnumbers)*2780)
grid on
grid minor
xlabel('Pseudo Rib Number','FontSize',15)
ylabel('Mass-Kg','FontSize',15)
legend('Rib Mass','Skin-Panel Mass','Combined Mass')
title('Mass of Pseudo rib and D-cell skin vs rib number','FontSize',15)
hold off
xlim([1 15])


stress=[];

for i=1:length(yDcell)
    stress(i)=q2calcdcell(i)/dcellpanelthicknesses(i,3) ;
end 

figure()
plot(yDcell,stress/(10^6))
hold on
yline(150,'LineWidth',2,'LineStyle','--')
title('Stress Variation Across D-Cell','FontSize',16)
legend('Calculated Stress','Tresca Shear Stress')
grid on
grid minor
xlabel('Semi-Span m','FontSize',16)
ylim([0 200])
ylabel('Stress-MPa','FontSize',16)


figure()
plot(yDcell,(dcellpanelthicknesses(1:end,3)*1000))
hold on
plot(yDcell,ceil((dcellpanelthicknesses(1:end,3)*1000)))
legend('Calculated','Actual')
title('D-Cell Panel Thickness Distribution','FontSize',16)
xlabel('Semi-Span m','FontSize',16)
ylabel('Thickness-mm','FontSize',16)
grid on
grid minor


dstresses=[];

for i=1:length(yDcell)
    dstresses(i)=q2calcdcell(i)/dcellpanelthicknesses(i,2);
end 

figure()
plot(yDcell,dstresses/(10^6))
hold on
yline(180)
hold off
ylim([0 200])




bucklingstressesDcell=[];

KS=[44.22,42.98,39.32,35.44,31.56,26.92,23.8,21,17.88,14.76];

for i=1:length(yDcell)
    bucklingstressesDcell(i)=KS(i)*E*(tdcell/arclengthDcell(HTchorddcell(i)))^2 ;
end 

figure()
plot(yDcell,bucklingstressesDcell/(10^6))

%-------------------------
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


%FUNCTION TO CALCULATE CG AT EACH STATION
function cg = calculate_cg(chordlength)
x=(0:1/1000 : chordlength);
summoments=zeros(1,length(x));
aerofoilcoordinates=zeros(1,length(x));


for i=1:length(x)
    aerofoilcoordinates(i)=chordlength*5*0.15*((0.2969*sqrt(x(i))) -(0.1260*x(i)) - (0.3516*(x(i))^2)    + (0.2843*x(i)^3)  - (0.1015*x(i)^4)  );
    summoments(i)=aerofoilcoordinates(i)*x(i) ;
end 

cg=sum(summoments)/(trapz(aerofoilcoordinates));
end 



%Function to calculate Front spar height at each station



%Function to calculate LH

function LH=calculate(V,n,W)




cmo_airfoil=-0.13 ;
ARwing=9;
Gamma_quarterchord=31.5;
Gamma_max=28.79;
CL_alpha_M_0=0.096*180/pi;
twist_degrees=-3;
Sexp_Sref=0.889;
F=1.294;
nairfoil=0.8;

M=V/299.517 ;
xw=30.86;
xcg=39.8394;
xh=74.26;
Sref=469.44;
SH=115.1504;
c=8.75;
lt=1.75; 
T=290970*4;

rho=(8.91/23.77)*1.225;

CL_alpha_M= (2*pi*ARwing*Sexp_Sref*F)/(  2+ sqrt(   4 + ((ARwing*ARwing*(1-(M*M)))/(nairfoil^2))*(1   +    ((tand(Gamma_max))^2)/(1-(M^2))      ) ));

%CL_alpha_M=9;

CM0w=(   (cmo_airfoil*(ARwing*ARwing*cosd(Gamma_quarterchord)/(ARwing+(2*cosd(Gamma_quarterchord)))))  - (0.01*twist_degrees)     )*(CL_alpha_M/CL_alpha_M_0) ;



% 
MO=CM0w*0.5*rho*V*V*Sref*c;

LH=(MO+((n*W*9.81)*(xcg-xw)) + (T*lt)   )/(xh-xcg);




for i=1:100
    LH= ( MO+((  (n*W*9.81) -LH   )*(xcg-xw)) + (T*lt)   ) /(xh-xcg);
end 


%method 2

%LH=(MO+((n*W*9.81)*(xcg-xw))+(T*lt))/(xh-xw) ;


end 




%----Function to calculate tfront and trear-----

function nonlinear = equations(vars,Torque,Fy,effectiveT,h1,h2,wbw,E)

tfront=vars(1);
trear=vars(2);


width=0.5*(h1-h2)/( sin(atan(   (h1-h2)/(2*wbw)      ))       ) ;

B1=(effectiveT*width/6)*(2 + (h2/h1)) + (tfront*h1/6);
B2=(effectiveT*width/6)*(2 + (h1/h2)) + (trear*h2/6);

Ix=(B1*h1*h1/2) + (B2*h2*h2/2);

nonlinear=[ (Fy*B1*h1/(Ix*tfront)) + (Torque/((h1+h2)*wbw*tfront)) - (8.1*E*( tfront^2/h1^2   )) ;
              (Fy*B2*h2/(Ix*trear)) - (Torque/((h1+h2)*wbw*trear)) - (8.1*E*( trear^2/h2^2   )) ];

end 



%-----function to calculate arc length---------


function arclength = arclengthDcell(chord)


NACA0015=[
0.0000     0.00000
0.0125     0.02367
0.0250     0.03268
0.0500     0.04443
0.0750     0.05250
0.1000     0.05853
0.1500     0.06682
0.2000     0.07172
0.2500     0.07427
0.3000     0.07502
0.4000     0.07254
0.5000     0.06617
0.6000     0.05704
0.7000     0.04580
0.8000     0.03279
0.9000     0.01810
0.9500     0.01008
1.0000     0.00158].*chord;




x=NACA0015(1:8,1);
y=NACA0015(1:8,2);


% Create the spline
pp = spline(x, y);

% Evaluate the spline at a fine set of points
t = linspace(x(1), x(end), 1000); % Fine grid for evaluation
spline_points = ppval(pp, t);

% Compute the derivative of the spline
spline_deriv = fnder(pp); % Derivative of the spline
deriv_points = ppval(spline_deriv, t);

% Compute the arc length using numerical integration
arclength = trapz(t, sqrt(1 + deriv_points.^2));

end 





%----Skin stringer Optimisation

function [panelareanostringers,panelareastringers,L,effectivethicknesst,thicknessnostringers,ribvolume,skinstringervolume,thicknessstringers,correctedstress,AS]=stringerskin(index,b,AS_btratioarray,ts_tratioarray,stressratioarray,HT_chord,bending_moment,E,y,NACA0015,wingboxheight,wingboxwidth)
D=[];
Tr=8e-3;  %initial rib thickness



farrar=0.75;

AS_btratio=AS_btratioarray(index);   %corresponding to farrar efficiency of 0.7
ts_tratio=ts_tratioarray(index);    %corresponding to farrar efficiency of 0.7
stressratio=stressratioarray(index);  %corresponding to above values of AS_bt and ts_t ratio


thicknessnostringers=[];
thicknessstringers=[];
N=[];
AS=[];
ts=[];
Nstringer=[];
L=[];
Lopt=[];

correctedstress=[];
effectivethicknesst=[];


panelareanostringers=[];
panelareastringers=[];

numberofstringers=[];

criticalshearbuckling=[];

for i=1:length(y)
    wingboxwidth(i)=0.5*HT_chord(i);
    wingboxheight(i)=min( HT_chord(i)*2*NACA0015(8,2)    , HT_chord(i)*2*NACA0015(14,2)        );
    N(i)=abs(bending_moment(i))/(wingboxheight(i)*wingboxwidth(i));
    thicknessnostringers(i)= (   (N(i) *wingboxwidth(i)^2)/(3.62*E)        )^(1/3) ;

    
    Nstringer(i)=abs(bending_moment(i))/(wingboxheight(i)*wingboxwidth(i));
    thicknessstringers(i)= (   (Nstringer(i) *b^2)/(3.62*E)        )^(1/3) ;

    AS(i)=AS_btratio*thicknessstringers(i)*b;
    ts(i)=ts_tratio*thicknessstringers(i);


    correctedstress(i)=Nstringer(i)*stressratio/thicknessstringers(i) ;


    L(i)=((farrar^2)*(Nstringer(i))*E)/(correctedstress(i))^2 ;

    D(i)=HT_chord(i)*2*NACA0015(9,2)  ;
    Lopt(i)=(4*(farrar^2)* (D(i))^2 *(Tr^2)*E  )/(Nstringer(i));

    panelareastringers(i)=(thicknessstringers(i) + (AS(i)/b))*b*(wingboxwidth(i)/b);
    panelareanostringers(i)=(wingboxwidth(i)*thicknessnostringers(i));

    numberofstringers(i)=(wingboxwidth(i)/b)-1 ;

    effectivethicknesst(i)=thicknessstringers(i) + (AS(i)/b);

    criticalshearbuckling(i)=8.1*E*(thicknessstringers(i)/b)^2 ;

end 



thicknessrib=5e-3;
Lchosen=L(length(L)/2);

y_Lweightarray=[0:Lchosen:y(length(y))];
x1=0;
y1=HT_chord(1);
x2=y(length(y));
y2=HT_chord(length(HT_chord));
x=y_Lweightarray;


    m = (y2 - y1) / (x2 - x1);
    c = y1 - m * x1;
    HTchord_Lweight = m * x + c;


Lweightarray=[];

for i=1:length(y_Lweightarray)
Lweightarray(i)=(0.5*HTchord_Lweight(i))*(HTchord_Lweight(i)*NACA0015(8,2)*2)*thicknessrib;
end 

ribvolume=sum(Lweightarray);
skinstringervolume=abs(trapz(panelareastringers,y));
%skinstringervolume=sum(panelareastringers)*0.01;



end 

function z = generate_step_inputs(x, y)
    % Round the values in y to the nearest integer
    y_rounded = round(y);
    
    % Initialize the z array with zeros
    z = zeros(size(y_rounded));
    
    % The first value in z is the same as the first value in y_rounded
    z(1) = y_rounded(1);
    
    % Iterate through the array and create step inputs
    for i = 2:length(y_rounded)
        if y_rounded(i) ~= y_rounded(i-1)
            z(i) = y_rounded(i);
        else
            z(i) = z(i-1);
        end
    end
end

%------------Function which transfers torque Load to D cell box

function nonlinear2=Torque2(variables,torque,Fy,frontsparheight,rearsparheight,wingboxwidth,G,E,effectivetfinal,Dcell_length,Dcellarea)

tfront=variables(1);
trear=variables(2);
q1=variables(3);
q2=variables(4);



Torque=abs(torque);
Shearforce=abs(Fy);
tdcell=5e-3;
wingboxarea=0.5*(frontsparheight+rearsparheight)*wingboxwidth ;


dthetadz_dcell=(1/(2*Dcellarea*G))*( q2*Dcell_length/tdcell    + (q2-q1-(Shearforce/2*frontsparheight))*(frontsparheight/tfront)          );
dthetadz_box=(1/(2*wingboxarea*G))* (  ((Shearforce/(2*frontsparheight))  + q1 - q2)*(frontsparheight/tfront)    +  (2*q1*wingboxwidth/effectivetfinal)  + (q1 - (Shearforce/(2*rearsparheight)))*(rearsparheight/trear)    )  ;


eq1=dthetadz_box-dthetadz_dcell ;
eq2=Torque-(2*wingboxarea*q1)-(2*Dcellarea*q2);
eq3=((q1 - (Shearforce/(2*rearsparheight)))/(trear)) - (8.1*E*((trear^2)/(rearsparheight^2))) ;
eq4=(q2-q1-(Shearforce/(2*frontsparheight)))/(tfront)  - (8.1*E*((tfront^2)/(frontsparheight^2))) ;

nonlinear2=[ eq1;
    eq2;
    eq3;
    eq4];


end 


function result = remove_negatives(arr)
    % Use logical indexing to filter out negative values
    result = arr(arr >= 0);
end