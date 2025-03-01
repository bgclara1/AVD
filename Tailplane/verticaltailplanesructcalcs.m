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


%----------------------------

engine_momentarm=20.5;
tailplane_momentarm=72-39 ;
maxthrust= 290000;
VTspan= 11.26 ;
VTarea=63.4;
taper=0.4 ;

VT_croot=(2*VTarea/((1+taper)*VTspan));
VT_ctip=taper*VT_croot ;

Maxload= maxthrust*engine_momentarm/tailplane_momentarm ;
plythickness=0.125e-03;
compressivestrength_ply=1.5e09;
shearstress_ply=8e07;




y=0:0.01:VTspan ;


VT_chord = trapezoidal(0,VT_croot,VTspan,VT_ctip,y);

VTlift=[];

for i=1:length(y)
    VTlift(i)=(4*Maxload/(pi*VTspan))*(1  - (y(i)^2)/(VTspan^2)    )^0.5 ;
end 

% figure()
% plot(y,VTlift,'r','LineWidth',2)
% hold on
% plot(y,-VTlift,'r','LineWidth',2)
% grid on
% grid minor
% xlabel('Semi Span m')
% ylabel('Lift N/m')
% title('Lift Distribution')



%----SHEAR and Bending Distribution--------

totaldist=VTlift;
yreversed=flip(y);
totaldistreversed=flip(totaldist);

shear_force_reversed = cumtrapz(yreversed,totaldistreversed);
bending_moment_reversed= cumtrapz(yreversed,shear_force_reversed);

shear_force = flip(shear_force_reversed);
bending_moment = flip(bending_moment_reversed);

 % figure()
 % plot(y,-shear_force,'r','LineWidth',2)
 % hold on
 % plot(y,shear_force,'r','LineWidth',2)
 % title('Shear Force')
 % grid on
 % grid minor
 % xlabel('Semi Span')
 % ylabel('Shear Force N/m')
 % hold off


 figure()
  plot(y,bending_moment,'r','LineWidth',2)
  hold on
  plot(y,-bending_moment,'r','LineWidth',2)
  title('Bending Moment')
  grid on
  grid minor
   xlabel('Semi Span')
 ylabel('Bending Moment-N')
   hold off


%---------Torque Load----------------

Torque=[];

for i=1:length(y)
    Torque(i)=VTlift(i)*(0.45-0.25)*VT_chord(i) ;
end 

    t_rev = flip(Torque);
    torque_rev = cumtrapz(flip(y),t_rev);
    T = flip(torque_rev);

figure()
plot(y,(-T),'r','LineWidth',2)
hold on
plot(y,T,'r','LineWidth',2)
ylabel('Torque-N')
xlabel('Semi-Span')
title('Torque Loading')
grid on
grid minor



%-----WingGeometry----------------

frontsparheight=[];
rearsparheight=[];
boxskinwidth=[];
boxheight=[];

NACA0012=[ 0.0000000 0.0000000
 0.0005839 0.0042603
 0.0023342 0.0084289
 0.0052468 0.0125011
 0.0093149 0.0164706
 0.0145291 0.0203300
 0.0208771 0.0240706
 0.0283441 0.0276827
 0.0369127 0.0311559
 0.0465628 0.0344792
 0.0572720 0.0376414
 0.0690152 0.0406310
 0.0817649 0.0434371
 0.0954915 0.0460489
 0.1101628 0.0484567
 0.1257446 0.0506513
 0.1422005 0.0526251
 0.1594921 0.0543715
 0.1775789 0.0558856
 0.1964187 0.0571640
 0.2159676 0.0582048
 0.2361799 0.0590081
 0.2570083 0.0595755
 0.2784042 0.0599102
 0.3003177 0.0600172
 0.3226976 0.0599028
 0.3454915 0.0595747
 0.3686463 0.0590419
 0.3921079 0.0583145
 0.4158215 0.0574033
 0.4397317 0.0563200
 0.4637826 0.0550769
 0.4879181 0.0536866
 0.5120819 0.0521620
 0.5362174 0.0505161
 0.5602683 0.0487619
 0.5841786 0.0469124
 0.6078921 0.0449802
 0.6313537 0.0429778
 0.6545085 0.0409174
 0.6773025 0.0388109
 0.6996823 0.0366700
 0.7215958 0.0345058
 0.7429917 0.0323294
 0.7638202 0.0301515
 0.7840324 0.0279828
 0.8035813 0.0258337
 0.8224211 0.0237142
 0.8405079 0.0216347
 0.8577995 0.0196051
 0.8742554 0.0176353
 0.8898372 0.0157351
 0.9045085 0.0139143
 0.9182351 0.0121823
 0.9309849 0.0105485
 0.9427280 0.0090217
 0.9534372 0.0076108
 0.9630873 0.0063238
 0.9716559 0.0051685
 0.9791229 0.0041519
 0.9854709 0.0032804
 0.9906850 0.0025595
 0.9947532 0.0019938
 0.9976658 0.0015870
 0.9994161 0.0013419
 1.0000000 0.0012600];

for i=1:length(y)
    frontsparheight(i)=2*0.058*VT_chord(i);
    rearsparheight(i)=2*0.036*VT_chord(i);
    boxskinwidth(i)=0.5*VT_chord(i);

    boxheight(i)=(frontsparheight(i)+rearsparheight(i))/2 ;
end 




%--------Number of Plies ---- Webs ----Shear ------45

frontwebpilenumber_45=[];
rearwebpilenumber_45=[];

shearforcefront=[];
shearforcerear=[];

for i=1:length(y)

    shearforcefront(i)= abs( abs(T(i))*1/(2*0.5*(frontsparheight(i)+rearsparheight(i))*boxskinwidth(i))    + abs(shear_force(i))/(2*frontsparheight(i))) ;
    shearforcerear(i)= abs(    abs(T(i))*1/(2*0.5*(frontsparheight(i)+rearsparheight(i))*boxskinwidth(i))    - abs(shear_force(i))/(2*rearsparheight(i)) ) ;
    
   
   frontwebpilenumber_45(i)=shearforcefront(i)/(plythickness*shearstress_ply);
   rearwebpilenumber_45(i)=shearforcerear(i)/( plythickness*shearstress_ply  );



end 


frontwebpilenumber_45calc=ceil(frontwebpilenumber_45);
frontwebpilenumber_45calc=makeOddEven(frontwebpilenumber_45calc);
frontwebpilenumber_45calc(length(frontwebpilenumber_45calc))=2;
frontwebpilenumber_45calc=2*frontwebpilenumber_45calc ;

rearwebpilenumber_45calc=ceil(rearwebpilenumber_45);
rearwebpilenumber_45calc=makeOddEven(rearwebpilenumber_45calc);
rearwebpilenumber_45calc(length(rearwebpilenumber_45calc))=2;
rearwebpilenumber_45calc=2*rearwebpilenumber_45calc ;



frontweb90_0=[];
rearweb90_0=[];


for i=1:length(y)
    frontweb90_0(i)=ceil(0.1*frontwebpilenumber_45calc(i))*2 ;
    

   rearweb90_0(i)=ceil(0.1*rearwebpilenumber_45calc(i))*2 ;


end 

figure()
plot(y,frontwebpilenumber_45calc)
hold on
plot(y,frontweb90_0)
title('Front Spar Ply Distribution','FontSize',16)
legend('-45/45 Ply Number','0/90 Ply Number')
ylabel('Ply Number','FontSize',16)
xlabel('Semi-Span m','FontSize',16)
grid on 
grid minor
hold off

figure()
plot(y,rearwebpilenumber_45calc)
hold on
plot(y,rearweb90_0)
title('Rear Spar Ply Distribution','FontSize',16)
legend('-45/45 Ply Number','0/90 Ply Number')
ylabel('Ply Number','FontSize',16)
xlabel('Semi-Span m','FontSize',16)
grid on 
grid minor
hold off

index=round(length(shearforcefront)/(7)) ;


VTfrontsparheight_sections=[frontsparheight(1),frontsparheight(1*index),frontsparheight(2*index),frontsparheight(3*index),frontsparheight(4*index),frontsparheight(5*index),frontsparheight(6*index)];
VTrearsparheight_sections=[rearsparheight(1),rearsparheight(1*index),rearsparheight(2*index),rearsparheight(3*index),rearsparheight(4*index),rearsparheight(5*index),rearsparheight(6*index)];

NXYfrontspar_sections=[shearforcefront(1),shearforcefront(1*index),shearforcefront(2*index),shearforcefront(3*index),shearforcefront(4*index),shearforcefront(5*index),shearforcefront(6*index)];
NXYrearspar_sections=[shearforcerear(1),shearforcerear(1*index),shearforcerear(2*index),shearforcerear(3*index),shearforcerear(4*index),shearforcerear(5*index),shearforcerear(6*index)];


%---------------Front Spar--------------------------------------------

%-------Section 1 front spar-------------


%theta=[-45,45,0,-45,45,0,-45,45,0,-45,45,90,-45,45,90,-45,45,90,45];
theta=transpose([-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
0
90
-45
45
90
-45
45
0
-45







]);
theta=makeMirrorArray(theta);
theta;

[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(1));

ratiosection1_frontspar=NXYfrontspar_sections(1)/Nxycrit ;




% %------Section 2 Front Spar--------------

%theta=[-45,45,0,-45,45,0,-45,45,0,-45,45,90,-45,45,90,-45,45,90,-45]
theta=transpose([-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45









]);

 
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(2));
 
 ratiosection2_frontspar=NXYfrontspar_sections(2)/Nxycrit ;





 %---Section 3 Front Spar --------------------------------------

% theta=[-45,45,0,-45,45,0,-45,45,90,-45,45,90,-45,45,0,-45];
theta=transpose([-45
45
0
90
-45
45
0
-45
45
90
-45
45
0
0
90
-45
45
0
-45
45
90
-45
45
0
-45









]);

 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(3));

 ratiosection3_frontspar=NXYfrontspar_sections(3)/Nxycrit ;



 %---Section 4 Front Spar --------------------------------------

%theta=[-45,45,90,-45,45,90,-45,45,0,-45,45,0,-45];
theta=transpose([-45
45
0
90
-45
45
0
0
90
-45
0
0
90
45
0
-45
45
90
-45
45
0
-45








]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(4));

 ratiosection4_frontspar=NXYfrontspar_sections(4)/Nxycrit;


 %---Section 5 Front Spar --------------------------------------

%theta=[-45,45,90,-45,45,90,-45,45,0,0,-45]
theta=transpose([-45
45
90
-45
45
90
45
45
0
90
0
-45
45
90
-45
45
0
-45





]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(5));

ratiosection5_frontspar=NXYfrontspar_sections(5)/Nxycrit;


 %---Section 6 Front Spar --------------------------------------

%theta=[-45,45,90,-45,45,0,-45]
theta=transpose([-45
45
90
-45
90
45
90
-45
0
0
-45
90
0
-45







]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(6));

ratiosection6_frontspar=NXYfrontspar_sections(6)/Nxycrit;




%----Section 7, front Spar-------------------------------


theta=transpose([-45
45
90
-45
90
0
-45
90
0
-45





]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections(7));

ratiosection7_frontspar=NXYfrontspar_sections(7)/Nxycrit;




%-------Checking the sections between layups---------

index2=round(length(shearforcefront)/(19)) ;


VTfrontsparheight_sections2=[];
Nxyfrontsparsections2=[];

VTfrontsparheight_sections2(1)=frontsparheight(1);
Nxyfrontsparsections2=shearforcefront(1);

for i=2:19
    VTfrontsparheight_sections2(i)=frontsparheight(index2*(i-1));
    Nxyfrontsparsections2(i)=shearforcefront(index2*(i-1));

end 


%----Front spar intersection checks

theta=transpose([-45
45
0
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45


]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTfrontsparheight_sections2(4));

ratio=Nxyfrontsparsections2(4)/Nxycrit ;

 
%-----------Rear Spar--------------------------------------------------

%------Rear Spar Section 1--------------

 
theta=transpose([-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
0
45
-45


]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(1));

ratiosection1_rearspar=NXYrearspar_sections(1)/Nxycrit;



%---------Rear Spar Section 2-------------------

theta=transpose([-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
0
-45
90
0
-45

]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(2));

ratiosection2_rearspar=NXYrearspar_sections(2)/Nxycrit;


%---------Rear Spar Section 3-------------------

theta=transpose([-45
45
0
0
90
-45
45
0
-45
45
90
-45
45
0
-45
45
90
-45




]);
theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(3));

ratiosection3_rearspar=NXYrearspar_sections(3)/Nxycrit;




%---------Rear Spar Section 4-------------------

theta=transpose([-45
45
0
90
-45
45
0
-45
45
90
-45
45
0
90
-45



]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(4));

ratiosection4_rearspar=NXYrearspar_sections(4)/Nxycrit;




%---------Rear Spar Section 5-------------------

theta=transpose([-45
45
0
90
-45
45
0
0
90
-45
-45
0
-45



]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(5));

ratiosection5_rearspar=NXYrearspar_sections(5)/Nxycrit;




%---------Rear Spar Section 6-------------------

theta=transpose([-45
45
90
-45
0
0
90
-45
0
-45


]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(6));

ratiosection6_rearspar=NXYrearspar_sections(6)/Nxycrit;


%---------Rear Spar Section 7-------------------

theta=transpose([-45
45
0
-45
45
90
-45



]);
 theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections(7));

ratiosection7_rearspar=NXYrearspar_sections(7)/Nxycrit;



%-------Checking the sections between layups---------

index2=round(length(shearforcerear)/(19)) ;


VTrearsparheight_sections2=[];
Nxyrearsparsections2=[];

VTrearsparheight_sections2(1)=rearsparheight(1);
Nxyrearsparsections2=shearforcerear(1);

for i=2:19
    VTrearsparheight_sections2(i)=rearsparheight(index2*(i-1));
    Nxyrearsparsections2(i)=shearforcerear(index2*(i-1));

end 


%----Rear spar intersection checks

theta=transpose([-45
-45
0
-45
45
90
-45
45
0
-45
45
90
-45
45
0
-45


]);

theta=makeMirrorArray(theta);
[Nxycrit]=ABD(theta,VTspan,VTrearsparheight_sections2(10));

ratio=Nxyrearsparsections2(10)/Nxycrit ;




%----Number of Plies ----Skin----Compression---0

numberofpileskin0=[];
numberofpilesksin9045=[];


for i=1:length(y)
    numberofpileskin0(i)=bending_moment(i)/(rearsparheight(i)*boxskinwidth(i)*plythickness*compressivestrength_ply) ;
    numberofpilesksin9045(i)=ceil(0.1*numberofpileskin0(i))*2 ;
end 

% 
 numberofpileskin0=ceil(numberofpileskin0);
 numberofpileskin0=makeOddEven(numberofpileskin0);
numberofpileskin0(length(numberofpileskin0))=2;
numberofpilesksin9045(length(numberofpilesksin9045))=2;
% 
 figure()
plot(linspace(0,VTspan,length(numberofpileskin0)),numberofpileskin0,'LineStyle','-')
hold on
plot(y,numberofpilesksin9045,'LineStyle','--')
legend('0 degrees Plies','+- 45 & 90 degrees Plies')
xlabel('Semi-Span m')
ylabel('Ply Number')
title('Skin Ply Number')
ylim([0 3])
grid on
grid minor
hold off





%-----Skin Sections------------------------

indexskin=round(length(bending_moment)/7) ;

wingboxskin_sections_width=[boxskinwidth(1),boxskinwidth(indexskin),boxskinwidth(2*indexskin),boxskinwidth(3*index),boxskinwidth(4*index),boxskinwidth(5*index),boxskinwidth(6*index)];

wingboxheightsections=[frontsparheight(1),frontsparheight(indexskin),frontsparheight(2*indexskin),frontsparheight(3*index),frontsparheight(4*index),frontsparheight(5*index),frontsparheight(6*index)];

bending_moment_sections=[bending_moment(1),bending_moment(indexskin),bending_moment(2*indexskin),bending_moment(3*index),bending_moment(4*index),bending_moment(5*index),bending_moment(6*index)];

Torqueskinsections=[T(1),T(indexskin),T(2*indexskin),T(3*index),T(4*index),T(5*index),T(6*index)];

VTsections=[y(1),y(indexskin),y(2*indexskin),y(3*indexskin),y(4*indexskin),y(5*indexskin),y(6*indexskin)];


b=0.3 ;
figure()
plot(y,boxskinwidth/b)


%-------Section 1----------------------

station=1;
b=0.3;
theta=transpose([-45
45
0
0
90
0
-45
45
0
0
90
0
0
-45
45
0
0
90
0
-45



]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX

%skinthickness=(plythickness*length(theta)*2);
%stringernumber=wingboxskin_sections_width(station)/b ;

%skinarea=skinthickness*wingboxskin_sections_width(station)
%stringerarea=0.5*b*skinthickness*wingboxskin_sections_width(station)/b 

%totalarea=skinarea+skinthickness


%---------Section 2---------------------------

station=2;
b=0.3;
theta=transpose([-45
45
0
0
90
-45
45
0
90
0
0
-45
45
0
0
90
0
-45


]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX




%---------Section 3---------------------------

station=3;
b=0.3;
theta=transpose([-45
45
0
0
90
-45
0
90
0
0
-45
45
0
90
0
-45

]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX


%---------Section 4---------------------------

station=4;
b=0.3;
theta=transpose([-45
45
0
0
-45
45
90
0
0
-45
45
90
0
-45

]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX



%---------Section 5---------------------------

station=5;
b=0.3;
theta=transpose([-45
45
0
0
-45
45
90
0
90
0
-45


]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX



%---------Section 6---------------------------

station=6;
b=0.3;
theta=transpose([-45
45
0
0
90
0
90
0
-45


]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX
% 
% 
% 
% %---------Section 7---------------------------
% 
station=7;
b=0.3;
theta=transpose([-45
45
90
0
-45


]);
theta=makeMirrorArray(theta);

NX=bending_moment_sections(station)/(frontsparheight(station)*wingboxskin_sections_width(station));
NXY=abs(Torqueskinsections(station))/(2*frontsparheight(station)*wingboxskin_sections_width(station));


[Nxcrit,Nxycrit]=ABD2(theta,VTspan,b);

ratiosection1_skin_NX=NX/Nxcrit
ratiosection1_skin_NXY=NXY/Nxycrit;
combinedloading=ratiosection1_skin_NXY^2 + ratiosection1_skin_NX



%---------Crushing Load Calculation----------------------
span_crushload=[0	0.4 	0.8 	1.2 	1.6	2 	2.4	2.8	3.2	3.605	4.01	4.415	4.82 5.2225	5.625	6.075	6.43 	6.8325	7.235	7.6375	8.04 	8.4425	8.845	9.2475	9.65	10.0525	10.445	10.8575	11.26 ];
thicknesscrushload=plythickness.*[40	38	38	38	36	34	34	34	32	30	30	30	28	26	24	24	22	20	20	20	18	16	14	12	10	10	10	10	10];
moment_crushload=linspace(bending_moment(1),bending_moment(length(bending_moment)),length(span_crushload));
AS_crushload=(0.5*0.3).*thicknesscrushload ;
boxheight_crushload=linspace(boxheight(1),boxheight(length(boxheight)),length(span_crushload));
chord_crushload=linspace(VT_chord(1),VT_chord(length(VT_chord)),length(span_crushload));
boxskinwidth_crushload=linspace(boxskinwidth(1),boxskinwidth(length(boxskinwidth)),length(span_crushload));

ribspacing=VTspan/6 ;
Ribpoints=0:ribspacing:VTspan;


%indexes based on spanpoints



riblocationspanindexes=[1,5,10,15,20,24,29];

riblocationspanindexes2=[1,162,322,483,644,805,967];

effectivethicknessrib=[];
moment_rib=[];
AS_rib=[];
boxheight_rib=[];
chord_rib=[];
span_rib=[0,1.61,3.2171,4.826,6.43,8.0429,9.65];
thicknessrib=plythickness.*[40,36,32,28,22,18,10];
boxskinwidth_rib=[];



for i=1:length(riblocationspanindexes2)
    % moment_rib(i)=moment_crushload(riblocationspanindexes(i));
    % AS_rib(i)=AS_crushload(riblocationspanindexes(i));
    % boxheight_rib(i)=boxheight_crushload(riblocationspanindexes(i));
    % chord_rib(i)=chord_crushload(riblocationspanindexes(i));
    % effectivethicknessrib(i)=thicknesscrushload(riblocationspanindexes(i))+(AS_rib(i)/0.3) ;
    % 
    % span_rib(i)=span_crushload(riblocationspanindexes(i)) ;
    % 
    % thicknessrib(i)=thicknesscrushload(riblocationspanindexes(i));
    % boxskinwidth_rib(i)=boxskinwidth_crushload(riblocationspanindexes(i));


   moment_rib(i)=bending_moment(riblocationspanindexes2(i));
   AS_rib(i)=0.5*0.3*thicknessrib(i);
   effectivethicknessrib(i)=thicknessrib(i)+(AS_rib(i)/0.3);
   boxheight_rib(i)=boxheight(riblocationspanindexes2(i));
   boxskinwidth_rib(i)=boxskinwidth(riblocationspanindexes2(i));
   chord_rib(i)=VT_chord(riblocationspanindexes2(i));


end 

thetarib_location_1=[-45 45 0 0 90 0 -45 45 0 0 90 0 0 -45 45 0 0 90 0 -45 ];
thetarib_location_1=makeMirrorArray(thetarib_location_1);


thetarib_location_2=[-45 45 0 0 90 -45 45 0 90 0 0 -45 45 0 0 90 0 -45];
thetarib_location_2=makeMirrorArray(thetarib_location_2);


thetarib_location_3=[-45 45 0 90 -45 0 90 0 0 -45 45 0 90 0 -45];
thetarib_location_3=makeMirrorArray(thetarib_location_3);



thetarib_location_4=[-45 45 0 0 -45 45 90 0 0 90 0 -45];
thetarib_location_4=makeMirrorArray(thetarib_location_4);


thetarib_location_5=[-45 45 0 -45 45 90 0 90 0 -45];
thetarib_location_5=makeMirrorArray(thetarib_location_5);



thetarib_location_6=[-45 45 0 90 0 -45];
thetarib_location_6=makeMirrorArray(thetarib_location_6);



thetarib_location_7=[-45 45 90 0 -45];
thetarib_location_7=makeMirrorArray(thetarib_location_7);


Eribarray=[ABD3(thetarib_location_1)/effectivethicknessrib(1),ABD3(thetarib_location_2)/effectivethicknessrib(2),ABD3(thetarib_location_3)/effectivethicknessrib(3),ABD3(thetarib_location_4)/effectivethicknessrib(4),ABD3(thetarib_location_5)/effectivethicknessrib(5),ABD3(thetarib_location_6)/effectivethicknessrib(6),ABD3(thetarib_location_7)/effectivethicknessrib(7)];


crushload_rib_locations=[];



for i=1:length(riblocationspanindexes)
    I=((chord_rib(i)*effectivethicknessrib(i)^3)/12)    + chord_rib(i)*effectivethicknessrib(i)*(boxheight_rib(i)/2)^2 ;
    
    crushload_rib_locations(i)=((moment_rib(i)^2)*ribspacing*boxheight_rib(i)*effectivethicknessrib(i)*chord_rib(i))/(2*Eribarray(i)*I*I);

end 

crushload_rib_locations2=[crushload_rib_locations,0];
span_rib2=[span_rib,11.26];


figure()
plot(span_rib2,crushload_rib_locations2,'LineWidth',2)
hold on
plot(span_rib2(1:end-1),crushload_rib_locations2(1:end-1),'LineStyle','none','Marker','x','MarkerSize',10)
grid on
grid minor
legend('Crush Load','Rib Location')
title('Crush Load Variation at Rib Locations','FontSize',12)
xlabel('Semi-Span m','FontSize',15)
ylabel('Crush Load N','FontSize',15)
xlim([0 11.26])


figure()
plot(span_rib,Eribarray)

kocoefficientarray=[0.026548673	0.03539823	0.115044248	0.194690265	0.265486726	0.283185841	0.300884956	0.300884956	0.309734513	0.309734513	0.318584071	0.327433628	0.336283186	0.336283186	0.336283186	0.345132743	0.345132743	0.353982301	0.371681416	0.371681416	0.380530973	0.380530973	0.389380531	0.398230088	0.398230088	0.415929204	0.433628319	0.442477876	0.451327434	0.451327434	0.469026549	0.469026549	0.442477876	0.495575221	0.495575221	0.495575221	0.495575221	0.522123894	0.548672566	0.592920354	0.610619469	0.637168142	0.699115044	0.752212389	0.814159292	0.920353982	1.061946903	1.159292035	1.274336283	1.274336283	1.389380531	1.389380531	1.442477876	1.522123894	1.628318584	1.716814159	1.716814159	1.805309735	1.911504425	2.053097345	2.150442478	2.283185841	2.389380531	2.513274336	2.628318584	2.725663717	2.831858407	2.946902655	4.769911504];
KOarray=[99.83416252	100.1658375	100.1658375	100.1658375	100.1658375	98.83913765	97.51243781	96.18573798	94.52736318	92.37147595	90.21558872	88.05970149	85.73797678	84.07960199	82.4212272	79.93366501	77.94361526	75.95356551	73.6318408	71.80762852	69.98341625	68.15920398	66.33499171	64.17910448	62.68656716	61.19402985	59.86733002	57.54560531	55.38971808	53.06799337	51.40961857	49.41956882	47.42951907	45.27363184	42.62023217	39.80099502	37.81094527	36.48424544	33.99668325	31.1774461	29.02155887	27.1973466	24.87562189	22.71973466	21.22719735	20.23217247	20.06633499	21.39303483	22.71973466	22.71973466	25.04145937	25.04145937	23.88059701	22.55389718	21.72470978	20.72968491	20.72968491	20.23217247	19.73466003	19.73466003	19.90049751	20.72968491	21.72470978	21.06135987	20.72968491	20.39800995	20.23217247	20.06633499	54.72636816];

%----Rib 1------------------------

station=1;
theta=[-45,45,0,0,90,0,0,-45,45,0,0,90,-45,45,0,0,90,0,0,-45,45,90,0,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib1=crushload_rib_locations(station)/Nxb 



%---------Rib 2-------------------------

station=2;
theta=[-45,45,0,0,90,0,0,-45,45,0,0,90,-45,45,0,90,0,-45,45,0,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib2=crushload_rib_locations(station)/Nxb 



%---------Rib 3-------------------------

station=3;
theta=[-45,0,45,0,90,0,0,-45,45,0,0,90,0,-45,45,90,0,0,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib3=crushload_rib_locations(station)/Nxb 




%--------Rib 4-----------------------------

station=4;
theta=[-45,0,45,0,90,0,0,-45,45,0,0,90,0,0,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib4=crushload_rib_locations(station)/Nxb 



%--------Rib 5------------------------------------


station=5;
theta=[-45,0,45,0,90,0,0,90,0,0,90,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib5=crushload_rib_locations(station)/Nxb 



%--------Rib 6------------------------------------


station=6;
theta=[-45,0,45,0,90,0,90,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib6=crushload_rib_locations(station)/Nxb



%-----------Rib 7-------------------------

station=7;
theta=[-45,45,0,90,-45];
theta=makeMirrorArray(theta);

[kocoefficient,D]=ABD4(0.5*chord_rib(station),boxheight_rib(station),theta);
bstation=boxskinwidth_rib(station);


calclosest=calculateclosestvalue(kocoefficient,kocoefficientarray,KOarray);
Nxb=(calclosest*sqrt(D(1,1)*D(2,2))/bstation^2) + (2*pi^2/bstation^2)*(D(1,2) + (2*D(3,3)));

stressratio_rib7=crushload_rib_locations(station)/Nxb


zerodegreepliesrib=[];

for i=1:7
    zerodegreepliesrib(i)=crushload_rib_locations(i)/(plythickness*compressivestrength_ply*boxheight_rib(i)*boxskinwidth_rib(i)) ;
end

figure()
plot(span_rib,ceil(zerodegreepliesrib))
hold on
plot(span_rib,ceil(0.1*zerodegreepliesrib*2))
hold on
plot(span_rib(1:end-1),ceil(zerodegreepliesrib(1:end-1)),'LineStyle','none','Marker','x','MarkerFaceColor','r','MarkerSize',16)
hold on
plot(span_rib(1:end-1),ceil(0.1*zerodegreepliesrib(1:end-1)*2),'LineStyle','none','Marker','x','MarkerFaceColor','r','MarkerSize',16)
legend('0 Degrees Plies','+/- 45 degrees Plies','Rib Locations','Rib Locations')
xlabel('Semi Span-m','FontSize',14)
ylabel('Number of Plies','FontSize',14)
title('Ribs Ply Distribution','FontSize',14)
ylim([0 2])
grid on
grid minor



%-------D-Cell--------------

tdcell=19e-3;
pseudoribnumbers=[1,3,5,7,9,11,13,15];
tresca=4.7e06;
E=70e09;

yDcell=linspace(0,VTspan,10);
VTchorddcell=trapezoidal(0,VT_croot,VTspan,VT_ctip,yDcell);

 a_bpsuedorib=zeros(length(VTchorddcell),length(pseudoribnumbers)) ;
 b_rt=zeros(length(VTchorddcell),length(pseudoribnumbers)) ;
 
 b_apseudorib=zeros(length(VTchorddcell),length(pseudoribnumbers)) ;
 a_rt=zeros(length(VTchorddcell),length(pseudoribnumbers)) ;



for j=1:length(pseudoribnumbers)
    a=VTspan/(pseudoribnumbers(j)+1);
for i=1:length(yDcell)
b=arclengthDcell(VTchorddcell(i));

R=VTchorddcell(i)*NACA0012(20,1);

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


KS=[58.1920904	54.6350365	27.37226277	18.39416058	13.02919708	10.83941606	9.525547445	7.99270073
53.78531073	57.26277372	28.79562044	19.48905109	14.01459854	11.27737226	9.854014599	8.321167883
49.71751412	60	30.87591241	20.69343066	14.7810219	11.93430657	10.18248175	8.540145985
46.10169492	59.7740113	34.59854015	22.77372263	15.76642336	12.48175182	10.51094891	8.97810219
42.14689266	57.40112994	37.66423358	24.52554745	16.53284672	13.35766423	11.05839416	9.635036496
38.98305085	51.52542373	40.40145985	26.27737226	17.95620438	14.12408759	11.60583942	9.854014599
35.14124294	46.66666667	44.23357664	28.90510949	19.9270073	15.32846715	12.48175182	10.2919708
31.29943503	38.8700565	45.42372881	31.75182482	21.89781022	16.64233577	13.35766423	11.16788321
28.36158192	33.89830508	40.33898305	36.02189781	25.40145985	18.17518248	14.45255474	12.26277372
25.53672316	29.03954802	35.14124294	38.75706215	29.45255474	22.00729927	17.73722628	14.23357664
];

 dcellpanelvolumes=[];
 dcellpanelthicknesses=[];
% 
 for j=1:length(pseudoribnumbers)
 for i=1:length(yDcell)
         dcellpanelthicknesses(i,j)=((arclengthDcell(VTchorddcell(i))^2)^(1/2))*(tresca/((KS(i,j)*E)))^(1/2);
 
         dcellpanelvolumes(i,j)=(arclengthDcell(VTchorddcell(i)))*dcellpanelthicknesses(i,j);
 end 
 end 
% 
 dcellpanelvolumes_differentribnumbers=[];
% 
 for i=1:length(pseudoribnumbers)
     dcellpanelvolumes_differentribnumbers(i)=abs(trapz(dcellpanelvolumes(1:end,i),yDcell));
 end 
% 
% 
 dcellribvolumes=[];
% 
 for j=1:length(pseudoribnumbers)
     ydcellribvolumes=linspace(0,VT_croot,pseudoribnumbers(j));
     VTchordRib=trapezoidal(0,VT_croot,11.26,VT_ctip,ydcellribvolumes);
     for i=1:length(ydcellribvolumes)
         dcellribvolumes(i,j)=arclengthDcell(VTchordRib(i))*(5e-3) ;
     end 
 
 end 
% 
 dcellvolumes_differentribnumbers=[];
% 
 for i=1:length(pseudoribnumbers)
 dcellvolumes_differentribnumbers(i)=sum(dcellribvolumes(1:end,i));
 end 
 
figure()
plot(pseudoribnumbers,dcellvolumes_differentribnumbers*2780)
hold on
plot(pseudoribnumbers,dcellpanelvolumes_differentribnumbers*2780)
hold on
plot(pseudoribnumbers,(dcellpanelvolumes_differentribnumbers+dcellvolumes_differentribnumbers)*2780)
grid on
grid minor
xlabel('Pseudo Rib Number','FontSize',15)
ylabel('Mass-Kg','FontSize',15)
legend('Rib Mass','Skin-Panel Mass','Combined Mass')
title('Mass of Pseudo rib and D-cell skin vs rib number','FontSize',15)
hold off
xticks(1:2:15);
xticklabels(1:2:15);
xlim([1 15])


 figure()
 plot(yDcell,(dcellpanelthicknesses(1:end,2)*1000))
 hold on
 plot(yDcell,ceil((dcellpanelthicknesses(1:end,2)*1000)))
 legend('Calculated','Actual')
 title('D-Cell Panel Thickness Distribution','FontSize',16)
 xlabel('Semi-Span m','FontSize',16)
 ylabel('Thickness-mm','FontSize',16)
 grid on
 grid minor



 %-----Weight Calculation-----VT

 %skin

skin_span=[0	0.4021	0.8043	1.2064	1.6086	2.0107	2.4129	2.815	3.2171	3.6193	4.0214	4.4236	4.8257	5.2279	5.63	6.0321	6.4343	6.8364	7.2386	7.6407	8.0429	8.445	8.8471	9.2493	9.6514	10.0536	10.4557	10.8579	11.26];
skin_thicknesspsan=(0.125e-3).*[40	38	38	38	36	34	34	34	32	30	30	30	28	26	24	24	22	20	20	20	18	16	14	12	10	10	10	10	10];
widthspan=linspace(boxskinwidth(1),boxskinwidth(length(boxskinwidth)),length(skin_span));
 
areaskin=[];

for i=1:length(skin_span)
    areaskin(i)=widthspan(i)*skin_thicknesspsan(i) ;
end 

skinweight=abs(trapz(skin_span,skin_thicknesspsan))*1550*2 


%frontsparweight

frontsparspan=linspace(0,11.26,29);
frontspar_thickness_span=(0.125e-3).*[60	60	60	58	56	56	54	52	50	48	48	46	44	42	40	38	36	34	32	30	28	26	24	22	20	20	20	20	20];
sparheightfront=linspace(frontsparheight(1),frontsparheight(length(frontsparheight)),length(frontsparheight)) ;

frontspararea=[];

for i=1:length(frontsparspan)
    frontspararea(i)=frontspar_thickness_span(i)*sparheightfront(i);
end 

frontsparweight=abs(trapz(frontsparspan,frontspararea))*1550




%Rear spar weight

rearspan=linspace(0,11.26,29);
rearspar_thickness_span=(0.125e-3).*[44	42	42	42	40	38	38	36	36	34	32	30	30	28	28	28	26	24	22	22	20	18	16	16	14	14	14	14	14];
sparheightrear=linspace(rearsparheight(1),rearsparheight(length(rearsparheight)),length(rearspan));

rearspararea=[];

for i=1:length(rearspan)
    rearspararea(i)=rearspar_thickness_span(i)*sparheightrear(i);
end


rearsparweight=abs(trapz(rearspan,rearspararea))*1550 


%ribweight

ribarea=[];

for i=1:length(span_rib)
    ribarea(i)=thicknessrib(i)*boxheight_rib(i)*0.5*chord_rib(i);
end 

ribweight=sum(ribarea)*1550


%---------------------------------------------------
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







%----ABD Matrix code---Cited from MATLAB

%A CODE
% Laminate definition (plies of equal thickness)

function Nxycrit = ABD(thetadt,span,b)

Nplies = length(thetadt);
E1   = 139.e9 ; % Pa
nu12 = 0.3 ;
E2   = 14.1e9  ; % Pa
G12  = 6.1e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
h_ply  = 0.125e-3;    % SI units, meters, ply thickness



a1 = -0.3e-6 ; % coefficients of thermal expansion
a2 = 28.1e-6 ;
deltaT = 1 ;





thetadb = fliplr(thetadt); % ply angles in degrees, from bottom          
h      = Nplies * h_ply ;
for i = 1:Nplies;
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end;
% Ply engineering properties (AS/3501)

% Q matrix (material coordinates)
 denom = 1 - nu12 * nu21 ;
Q11 = E1 / denom        ;
Q12 = nu12 * E2 / denom ;
Q22 = E2 / denom        ;
Q66 = G12               ;
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;
%Q=[20 .7 0; .7 2 0; 0 0 .7]
a = [a1 a2 0]' ;
S=inv(Q) ;
% Qbar matrices (laminate coordinates) and contributions to
% ABD matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % ply i angle in radians, from bottom
  m = cos(theta) ;
  n = sin(theta) ;
  T = [ m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)];
  Qbar = inv(T) * Q * (inv(T))' ;
  Sbar = inv(Qbar);

   abar = T' * a ;

  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);

   NT = NT + Qbar * abar * h_ply * deltaT ;
  MT = MT + Qbar * abar * h_ply * zbar(i) * deltaT ; 
end;
Qbar;
Sbar;
A;
B;
D;
ABD    = [A B; B D] ;

D11=D(1,1);
D12=D(1,2);
D22=D(2,2);
D33=D(3,3);

DOfactor=(D12+(2*D33))/(D11*D22)^(0.5);
D22_D11factor=(span/b)*((D22/D11)^(0.25));

D11_D22factor=(b/span)*((D11/D22)^(0.25));
phifactor=sqrt(D11*D22)/D33 ;

K=(2*D33 + D12)/sqrt(D11*D22);
Nxycrit=(4/(b^2))*((D11*(D22^3))^(1/4))*(8.125+(5.045*K));



end 





function [Nxcrit,Nxycrit] = ABD2(thetadt,span,b)

Nplies = length(thetadt);
E1   = 139.e9 ; % Pa
nu12 = 0.3 ;
E2   = 14.1e9  ; % Pa
G12  = 6.1e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
h_ply  = 0.125e-3;    % SI units, meters, ply thickness



a1 = -0.3e-6 ; % coefficients of thermal expansion
a2 = 28.1e-6 ;
deltaT = 1 ;





thetadb = fliplr(thetadt); % ply angles in degrees, from bottom          
h      = Nplies * h_ply ;
for i = 1:Nplies;
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end;
% Ply engineering properties (AS/3501)

% Q matrix (material coordinates)
 denom = 1 - nu12 * nu21 ;
Q11 = E1 / denom        ;
Q12 = nu12 * E2 / denom ;
Q22 = E2 / denom        ;
Q66 = G12               ;
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;
%Q=[20 .7 0; .7 2 0; 0 0 .7]
a = [a1 a2 0]' ;
S=inv(Q) ;
% Qbar matrices (laminate coordinates) and contributions to
% ABD matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % ply i angle in radians, from bottom
  m = cos(theta) ;
  n = sin(theta) ;
  T = [ m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)];
  Qbar = inv(T) * Q * (inv(T))' ;
  Sbar = inv(Qbar);

   abar = T' * a ;

  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);

   NT = NT + Qbar * abar * h_ply * deltaT ;
  MT = MT + Qbar * abar * h_ply * zbar(i) * deltaT ; 
end;
Qbar;
Sbar;
A;
B;
D;
ABD    = [A B; B D] ;

D11=D(1,1);
D12=D(1,2);
D22=D(2,2);
D33=D(3,3);


Nxcrit=2*((pi/b)^2)*(    sqrt(D11*D22)  + D12 + (2*D33));

DOfactor=(D12+(2*D33))/(D11*D22)^(0.5);
D22_D11factor=(span/b)*((D22/D11)^(0.25));

D11_D22factor=(b/span)*((D11/D22)^(0.25));
phifactor=sqrt(D11*D22)/D33 ;

K=(2*D33 + D12)/sqrt(D11*D22);
Nxycrit=(4/(b^2))*((D11*(D22^3))^(1/4))*(8.125+(5.045*K));


end 




function [Ecalcwithouth] = ABD3(thetadt)

Nplies = length(thetadt);
E1   = 139.e9 ; % Pa
nu12 = 0.3 ;
E2   = 14.1e9  ; % Pa
G12  = 6.1e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
h_ply  = 0.125e-3;    % SI units, meters, ply thickness



a1 = -0.3e-6 ; % coefficients of thermal expansion
a2 = 28.1e-6 ;
deltaT = 1 ;





thetadb = fliplr(thetadt); % ply angles in degrees, from bottom          
h      = Nplies * h_ply ;
for i = 1:Nplies;
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end;
% Ply engineering properties (AS/3501)

% Q matrix (material coordinates)
 denom = 1 - nu12 * nu21 ;
Q11 = E1 / denom        ;
Q12 = nu12 * E2 / denom ;
Q22 = E2 / denom        ;
Q66 = G12               ;
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;
%Q=[20 .7 0; .7 2 0; 0 0 .7]
a = [a1 a2 0]' ;
S=inv(Q) ;
% Qbar matrices (laminate coordinates) and contributions to
% ABD matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % ply i angle in radians, from bottom
  m = cos(theta) ;
  n = sin(theta) ;
  T = [ m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)];
  Qbar = inv(T) * Q * (inv(T))' ;
  Sbar = inv(Qbar);

   abar = T' * a ;

  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);

   NT = NT + Qbar * abar * h_ply * deltaT ;
  MT = MT + Qbar * abar * h_ply * zbar(i) * deltaT ; 
end;
Qbar;
Sbar;
A;
B;
D;
ABD    = [A B; B D] ;

D11=D(1,1);
D12=D(1,2);
D22=D(2,2);
D33=D(3,3);

term1=A(1,1);
term2=A(1,2)*(  A(2,3)*A(1,3)  - A(1,2)*A(3,3)       )/( A(2,2)*A(3,3)  - A(2,3)^2   ) ;
term3=A(1,3)*(   (-A(1,3)/(A(3,3)))    +      (A(2,3)*A(1,2)*A(3,3)     -    A(2,3)*A(2,3)*A(1,3)  )/(A(2,2)*A(3,3)*A(3,3) -   A(2,3)*A(2,3)*A(3,3)        )       ) ;

Ecalcwithouth=term1+term2+term3 ;

Ecalcwithouth=A(1,1)-(A(1,2)^2)/(A(2,2))


end 






function [Kocoefficient,D] = ABD4(wingboxwidth,wingboxheight,thetadt)

Nplies = length(thetadt);
E1   = 139.e9 ; % Pa
nu12 = 0.3 ;
E2   = 14.1e9  ; % Pa
G12  = 6.1e9  ; % Pa
nu21 = nu12 * E2 / E1 ;
h_ply  = 0.125e-3;    % SI units, meters, ply thickness



a1 = -0.3e-6 ; % coefficients of thermal expansion
a2 = 28.1e-6 ;
deltaT = 1 ;





thetadb = fliplr(thetadt); % ply angles in degrees, from bottom          
h      = Nplies * h_ply ;
for i = 1:Nplies;
  zbar(i) = - (h + h_ply)/2 + i*h_ply;
end;
% Ply engineering properties (AS/3501)

% Q matrix (material coordinates)
 denom = 1 - nu12 * nu21 ;
Q11 = E1 / denom        ;
Q12 = nu12 * E2 / denom ;
Q22 = E2 / denom        ;
Q66 = G12               ;
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66] ;
%Q=[20 .7 0; .7 2 0; 0 0 .7]
a = [a1 a2 0]' ;
S=inv(Q) ;
% Qbar matrices (laminate coordinates) and contributions to
% ABD matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
NT = zeros(3,1);
MT = zeros(3,1);
for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % ply i angle in radians, from bottom
  m = cos(theta) ;
  n = sin(theta) ;
  T = [ m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n (m^2 - n^2)];
  Qbar = inv(T) * Q * (inv(T))' ;
  Sbar = inv(Qbar);

   abar = T' * a ;

  A = A + Qbar * h_ply;
  B = B + Qbar * h_ply * zbar(i); 
  D = D + Qbar * (h_ply * zbar(i)^2  + h_ply^3 / 12);

   NT = NT + Qbar * abar * h_ply * deltaT ;
  MT = MT + Qbar * abar * h_ply * zbar(i) * deltaT ; 
end;
Qbar;
Sbar;
A;
B;
D;
ABD    = [A B; B D] ;

D11=D(1,1);
D12=D(1,2);
D22=D(2,2);
D33=D(3,3);

Kocoefficient=(wingboxheight/wingboxwidth)*(D22/D11)^(1/4) ;




end 








function outputArray = makeOddEven(inputArray)
    % makeOddEven: Converts all odd numbers in the input array to even numbers.
    % Input: inputArray - An array of integers
    % Output: outputArray - An array where all odd numbers are made even
    
    % Initialize the output array with the same size as the input array
    outputArray = inputArray;
    
    % Loop through each element in the array
    for i = 1:numel(inputArray)
        % Check if the current element is odd
        if mod(inputArray(i), 2) == 1
            % If odd, add 1 to make it even
            outputArray(i) = inputArray(i) + 1;
        end
    end
end


function mirroredArray = makeMirrorArray(inputArray)
    % makeMirrorArray: Creates a mirror image of the input array.
    % Input: inputArray - An array of numbers
    % Output: mirroredArray - The mirror image of the input array
    
    % Reverse the input array
    reversedArray = fliplr(inputArray);

    reversedArray(1)=-inputArray(length(inputArray));
    
    % Concatenate the original array with the reversed array
    mirroredArray = [inputArray, reversedArray];
end








function calclosest = calculateclosestvalue(a,calccoefficientarray,KOarray)

error=[];

for i=1:length(calccoefficientarray)
    error(i)=abs(calccoefficientarray(i)-a);
end 


calclosest=KOarray(find(error==min(error))) ;
calclosest=calclosest(1);


end 




function arclength = arclengthDcell(chord)


NACA0012=[0.0000000 0.0000000
 0.0005839 0.0042603
 0.0023342 0.0084289
 0.0052468 0.0125011
 0.0093149 0.0164706
 0.0145291 0.0203300
 0.0208771 0.0240706
 0.0283441 0.0276827
 0.0369127 0.0311559
 0.0465628 0.0344792
 0.0572720 0.0376414
 0.0690152 0.0406310
 0.0817649 0.0434371
 0.0954915 0.0460489
 0.1101628 0.0484567
 0.1257446 0.0506513
 0.1422005 0.0526251
 0.1594921 0.0543715
 0.1775789 0.0558856
 0.1964187 0.0571640
 0.2159676 0.0582048
 0.2361799 0.0590081
 0.2570083 0.0595755
 0.2784042 0.0599102
 0.3003177 0.0600172
 0.3226976 0.0599028
 0.3454915 0.0595747
 0.3686463 0.0590419
 0.3921079 0.0583145
 0.4158215 0.0574033
 0.4397317 0.0563200
 0.4637826 0.0550769
 0.4879181 0.0536866
 0.5120819 0.0521620
 0.5362174 0.0505161
 0.5602683 0.0487619
 0.5841786 0.0469124
 0.6078921 0.0449802
 0.6313537 0.0429778
 0.6545085 0.0409174
 0.6773025 0.0388109
 0.6996823 0.0366700
 0.7215958 0.0345058
 0.7429917 0.0323294
 0.7638202 0.0301515
 0.7840324 0.0279828
 0.8035813 0.0258337
 0.8224211 0.0237142
 0.8405079 0.0216347
 0.8577995 0.0196051
 0.8742554 0.0176353
 0.8898372 0.0157351
 0.9045085 0.0139143
 0.9182351 0.0121823
 0.9309849 0.0105485
 0.9427280 0.0090217
 0.9534372 0.0076108
 0.9630873 0.0063238
 0.9716559 0.0051685
 0.9791229 0.0041519
 0.9854709 0.0032804
 0.9906850 0.0025595
 0.9947532 0.0019938
 0.9976658 0.0015870
 0.9994161 0.0013419
 1.0000000 0.0012600].*chord;




x=NACA0012(1:20,1);
y=NACA0012(1:20,2);


% Create the spline
pp = spline(x, y);

% Evaluate the spline at a fine set of points
t = linspace(x(1), x(end), 1000); % Fine grid for evaluation
spline_points = ppval(pp, t);

% Compute the derivative of the spline
spline_deriv = fnder(pp); % Derivative of the spline
deriv_points = ppval(spline_deriv, t);

% Compute the arc length using numerical integration
arclength = 2*trapz(t, sqrt(1 + deriv_points.^2));

end 
