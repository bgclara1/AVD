close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            SET DEFAULT PLOTTING PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',15)
set(groot,'defaulttextfontsize',15)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           INERTIAL LOAD DISCRETISATION AND PLOTTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lTot = 80;
DiscInterval = 5;
fuselageDisc = 0:DiscInterval:lTot;
cabinLow = 6.23;
cabinHigh = 60.48;
numberOfPoints = lTot/DiscInterval;

fuselagePoints = [];
fuselageWeight = 29744*9.81;
for i = 1:lTot
    fuselagePoints(i) = fuselageWeight/lTot;
    xDiscr(i) = i-1;
end

payloadPoints = [];
payloadWeight = 51879*9.81;

for i = 1:lTot
    if i<=6
        payloadPoints(i) = 0;
    elseif i>=61 
        payloadPoints(i) = 0;
    else 
        payloadPoints(i) = payloadWeight/lTot;
    end
end


furnishingPoints = [];
furnishingWeight = 17227*9.81;

for i = 1:lTot
    if i<=6
        furnishingPoints(i) = 0;
    elseif i>=61 
        furnishingPoints(i) = 0;
    else 
        furnishingPoints(i) = furnishingWeight/lTot;
    end
end

fuelPoints = [];
fuelWeight= 21000*9.81;

for i = 1:lTot
    if i<=31
        fuelPoints(i) = 0;
    elseif i>=39 
        fuelPoints(i) = 0;
    else 
        fuelPoints(i) = fuelWeight/lTot;
    end
end

components = {
    "Horizontal Tailplane", 2998, 75;
    "Vertical Tailplane", 1552.8, 72;
    "Nose Landing Gear", 486.52, 2;
    "Fuel Systems", 936, 46;
    "Flight Controls", 927, 47;
    "Installed APU", 1496, 71;
    "Instruments", 796.77, 3;
    "Hydraulic System", 271.9, 46;
    "Electronic System", 1037, 37;
    "Avionics", 697.64, 3;
    "Air Conditioning", 3593, 33; 
    "Anti Icing System", 757.9, 42;
    "Handling Gear", 113.3685, 33;
};

for i = 1:size(components)
    weights(i) = components{i, 2} * 9.81; % Mass (kg) * 9.81 to get weight (N)
    positions(i) = components{i, 3};      % x_cg position (m)
end

[uniquePositions, ~, idx] = unique(positions); % Unique positions and their indices
aggregatedWeights = accumarray(idx, weights); % Sum weights for each unique position



total = [];
for i = 1:lTot
    total(i) = fuselagePoints(i) + payloadPoints(i) + furnishingPoints(i) + fuelPoints(i);
end

%for i = 1:lTot
 %   total(i) = fuselagePoints(i)
%end



total(75) = total(75) + 2998*9.81; %HT
total(72) = total(72) + 1552.8*9.81; %VT
total(2) = total(2) + 486.52*9.81; %nose landing gear
total(46) = total(46) + 936*9.81; %fuel sus
total(47) = total(47) + 927*9.81; %flight controls
total(3) = total(3) + 796.77*9.81; %instruments
total(71) = total(71) + 1496*9.81; %APU
total(37) = total(37) + 1037*9.81; %Elec 
total(42) = total(42)+ 757.9*9.81; %anti icing


figure;
plot(xDiscr, total*-1, 'b-', 'LineWidth', 1.5);

xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Load (N/m)','FontWeight','bold');
title('Inertial Loads');
grid minor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           AIR LOAD PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lTot = 79;
Sref = 115.1508; %HT
SrefWing = 469.44;
frontSparWP = 0.2;
backSparWP = 0.8;
frontSparTP = 0.2;
chordHT = 8.25;

wingWidth = 13.6; %chord
WingLE = 28.4779;
xPosFSW = WingLE+0.2*wingWidth;
xPosBSW = WingLE + 0.8*wingWidth;
xPosFST = 74.86; % Nadir

Vd = 192.4;
Va = 92.4;

rho = 1.225;
rhoC = 0.4592; % C cruise

MAC = 8.75;
Cm = -0.14;

%------------------- 3.75 --------------------------------
MoW = 3.75*Cm*(0.5*rhoC*Vd^2*SrefWing*MAC);
LHT = (-330800*9.81*(38.8831-30.86)+MoW)/(71.1-30.86)

A = [ xPosFSW	xPosBSW ;
        1	1];
B = [xPosFST*LHT*-1; LHT*-1];


X = linsolve(A,B);
forceRF = X(1)
forceRR = X(2)



aero =[];
aerofuselageDisc = fuselageDisc;
aerofuselageDisc(end+1) = xPosBSW;
aerofuselageDisc(end+1) = xPosFSW;
aerofuselageDisc(end+1) = xPosFST;

aerofuselageDisc = sort(aerofuselageDisc);

airXDiscr = xDiscr;
aero = zeros(lTot+1);
aero(39) = forceRR; % back wing spar
aero(31) = forceRF; %front wing spar
aero(71) = LHT; % HT ac
aero = aero(:,1);





%-----------------------------------------------------

%------------------- - 1.5 --------------------------------

MoW15 = -1.5*Cm*(0.5*rhoC*Va^2*SrefWing*MAC);
LHT15 = (-330800*9.81*(38.8831-30.86)+MoW15)/(71.1-30.86);

A = [ xPosFSW	xPosBSW ;
        1	1];
B = [xPosFST*LHT15*-1; LHT15*-1];

X15 = linsolve(A,B);
forceRF15 = X15(1);
forceRR15 = X15(2);

aero15 = zeros(lTot+1);
aero15(39) = -forceRR15; % back wing spar
aero15(31) = -forceRF15; %front wing spar
aero15(71) = -LHT15; % HT ac
aero15 = aero15(:,1);


%-----------------------------------------------------

%{
figure;
plot(xDiscr, total*-1, 'b-', 'LineWidth', 1.5);

xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Load (N/m)','FontWeight','bold');
title('Inertial Loads');
grid minor

hold on;
stem(xDiscr,aero, 'LineWidth', 3, 'Marker', 'none', 'Color', 'r');
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Load (N/m)','FontWeight','bold');

title('Total Loads');

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SHEAR FORCE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



totalInertialLoad = sum(total);

for i = 1:length(total)
    inertialMomentArray(i) = total(i)*(i);
end

inertialMoment = sum(inertialMomentArray);

A = [ 31	39 ;
        1	1];
B = [inertialMoment; totalInertialLoad];


X = linsolve(A,B);
forceRFInert = X(1);
forceRRInert = X(2);

MomentRFInert = forceRFInert*xPosFSW;
MomentRRInert = forceRRInert*xPosBSW;

SparReactInertLoadTot =  total*-1;
SparReactInertLoadTot(31) = SparReactInertLoadTot(31) - forceRFInert;
SparReactInertLoadTot(39) = SparReactInertLoadTot(39) + forceRRInert;

SFInertial = [];
SFInertial(1) = SparReactInertLoadTot(1);
for i = 2:length(SparReactInertLoadTot)
    SFInertial(i) = SFInertial(i-1) + SparReactInertLoadTot(i);
end

SFAero = [];
SFAero(1) = aero(1);
for i = 2:length(aero)
    SFAero(i) = SFAero(i-1) + aero(i);
end

SF375 = [];
SF375(1) = SFAero(1) + SFInertial(1);
for i = 2:length(SFAero)
    SF375(i) = SFAero(i) + SFInertial(i);
end

SFaero15 = [];
SFaero15(1) = aero15(1);
for i = 2:length(aero15)
    SFaero15(i) = SFaero15(i-1) + aero15(i);
end

SF15 = [];
SF15(1) = SFaero15(1) + SFInertial(1);
for i = 2:length(SFaero15)
    SF15(i) = SFaero15(i) + SFInertial(i);
end


%============= OEI =========================

RC = 8.05; 
TC = 0.3 * RC;
MAC = (2/3) * ((RC + TC) - ((RC * TC) / (RC + TC)))

RC = 8.05;
TC = 0.3 * RC; 
Span = 11.26; 

MAC = (2/3) * ((RC + TC) - ((RC * TC) / (RC + TC)));
d_MAC = (Span / 3) * ((RC + 2 * TC) / (RC + TC));
Diam = 6.485;
wing_offset = 1.36;
d_wing_MAC = d_MAC + Diam/2 + wing_offset;

OEI_Force = 1.8979*1e+5;

OEI_Moment = d_wing_MAC * OEI_Force;
SemiSpan = 32.5;

OEI_Aileron = OEI_Moment/(SemiSpan/2);



%OEI_Moment = OEI_Force*OEI_Arm;

A = [xPosFSW xPosBSW; 
     1  1]; 
B = [OEI_Moment; OEI_Aileron];

X = linsolve(A, B);
forceRFInertOEI = X(1);
forceRRInertOEI = X(2);

SparReactInertLoadTot = total * -1; 

SparReactInertLoadTot(31) = SparReactInertLoadTot(31) - forceRFInertOEI;
SparReactInertLoadTot(39) = SparReactInertLoadTot(39) - forceRRInertOEI;

%SparReactInertLoadTot(OEI_Position) = SparReactInertLoadTot(OEI_Position) + OEI_Thrust;

SF_OEI = [];
SF_OEI(1) = SparReactInertLoadTot(1);
for i = 2:length(SparReactInertLoadTot)
    SF_OEI(i) = SF_OEI(i-1) + SparReactInertLoadTot(i);
end

%---------------- Landing --------------------------
    
MaxLandingMassFrac = 0.85;
MoW15 = 2.18*Cm*(0.5*rhoC*Va^2*SrefWing*MAC);
LHTland = (MaxLandingMassFrac*-330800*9.81*(38.8831-30.86)+MoW15)/(71.1-30.86);

A = [ xPosFSW	xPosBSW ;
        1	1];
B = [xPosFST*LHT15*-1; LHT15*-1];

Xland = linsolve(A,B);
forceRFland = Xland(1);
forceRRland = Xland(2);

aeroland = zeros(lTot+1);
aeroland(39) = forceRR15; % back wing spar
aeroland(31) = forceRF15; %front wing spar
aeroland(71) = LHTland; % HT ac
aeroland = aeroland(:,1);



SFaeroland = [];
SFaeroland(1) = aeroland(1);
for i = 2:length(aeroland)
    SFaeroland(i) = SFaeroland(i-1) + aeroland(i);
end

SFland = [];
SFland(1) = SFaeroland(1) + SFInertial(1);
for i = 2:length(SFaeroland)
    SFland(i) = SFaeroland(i) + SFInertial(i);
end



figure;
plot(xDiscr, SFInertial, 'LineWidth', 1.5);
hold on;
plot(xDiscr, SFaero15, 'LineWidth', 1.5);
hold on;
plot(xDiscr, SFAero, 'LineWidth', 1.5);
hold on;
plot(xDiscr, SF_OEI, 'r-', 'LineWidth', 1.5); % OEI case
hold on;
plot(xDiscr, SF375, 'LineWidth', 1.5)
hold on;
plot(xDiscr, SF15, 'LineWidth', 1.5)
hold on;
plot(xDiscr, SFland, 'LineWidth', 1.5)
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Shear Force (N/m)','FontWeight','bold');
title('Shear Force along Fuselage');
legend('Inertial', 'Aero -1.5', 'aero 3.75','OEI Case', '3.75', '-1.5', 'landing');
%legend('Inertial',  '3.75', '-1.5', 'landing');
grid minor

%------------ both SF --------------------

SFTot375 = SFAero+SFInertial;
SFTot15 = SFaero15+SFInertial;
%{
figure;
plot(xDiscr, SFTot375, 'LineWidth', 1.5);
hold on;
plot(xDiscr, SFTot15, 'LineWidth', 1.5);
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Shear Force (N/m)','FontWeight','bold');
title('Shear Force along Fuselage');
grid minor

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          INERTIAL BENDING MOMENT PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dBM = [];
dBM(1) = 0;
for i = 2:length(SFInertial)
    dBM(i) = ((SFInertial(i-1)*(xDiscr(i)-xDiscr(i-1))+SFInertial(i)*(xDiscr(i)-xDiscr(i-1)))/2);
end

BM = [];
BM(1) = dBM(1);
for i = 2:length(dBM)
    BM(i) = BM(i-1) + dBM(i);
end

BMInert = cumtrapz(SFInertial);

figure;
plot(xDiscr, BM, 'b-', 'LineWidth', 1.5);
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Bending Moment(Nm)','FontWeight','bold');
title('Bending Moment along Fuselage');
grid minor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          AIR BENDING MOMENT PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BMAero = cumtrapz(SFAero);
%{
figure;
plot(xDiscr, BMAero, 'b-', 'LineWidth', 1.5);
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Bending Moment(Nm)','FontWeight','bold');
title('Bending Moment along Fuselage');
grid minor
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          AIR AND INERTIAL BENDING MOMENT PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BMTot = BMAero+BMInert;
%{
plot(xDiscr, BMTot, 'b-', 'LineWidth', 1.5);
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Bending Moment(Nm)','FontWeight','bold');
title('Bending Moment along Fuselage');
grid minor
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           GEAR LOAD PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



MTOW = 330800;
MaxLandingMass = MTOW*0.85;
landingLoad = 2.18*MaxLandingMass*9.81;
GearReaction = landingLoad;

gearXDiscr = xDiscr;
gear = zeros(lTot);
gear(39) = GearReaction; % back wing spar

  
%{
figure;
plot(gearXDiscr, gear, 'b-', 'LineWidth', 1.5);
xlabel('Fuselage Position (m)','FontWeight','bold');
ylabel('Landing Gear Load (N/m)','FontWeight','bold');
title('Gear Reaction Load');
grid minor

%}






