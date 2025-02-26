function [total, comboSF] = combineOEI(xDiscr,InertialLoads,LHTOEI)

    % all inertial and LHT loads are like before but we're gonna add a load
    % at midspan to act as the reaction force to the rolling moment from
    % the rudder restorative force. then work out the spar reaction forces
    % from the total loads incl the aileron force.

    Load = InertialLoads*1;


    Load(71) = Load(71)+ LHTOEI/2;
    Load(75) = Load(75)+ LHTOEI/2;

    RC = 8.05; 
    TC = 0.3 * RC;
    MAC = (2/3) * ((RC + TC) - ((RC * TC) / (RC + TC)));
    
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
    
    OEI_Aileron = OEI_Moment/(SemiSpan/2)

    Load(35) = Load(35) + OEI_Aileron;

    totalLoad = sum(Load);

    for i = 1:length(xDiscr)
        inertialMoments(i) = Load(i) * xDiscr(i);
    end
    

    totalMoment = sum(inertialMoments);

    A = [ 31	39 ;
            1	1];
    B = [totalMoment; totalLoad];
    
    X = linsolve(A,B);
    forceRF = X(1)*-1;
    forceRR = X(2)*-1;


    sparReaction = zeros(1,80);
    sparReaction(31) = forceRF;
    sparReaction(39) = forceRR;

    for i = 1:length(xDiscr)
        total(i) = Load(i)  + sparReaction(i);
    end
    

    comboSF = total(1);
    for i = (2:length(total))
        comboSF(i) = comboSF(i-1) + total(i);
    end



end