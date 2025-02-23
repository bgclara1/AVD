function SF_OEI = getOEI(xDiscr, total) 
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
    
    OEI_Aileron = OEI_Moment/(SemiSpan/2);

    
    wingWidth = 13.6; %chord
    WingLE = 28.4779;
    xPosFSW = WingLE+0.2*wingWidth;
    xPosBSW = WingLE + 0.8*wingWidth;
    
    
    %OEI_Moment = OEI_Force*OEI_Arm;
    
    A = [xPosFSW xPosBSW; 
         1  1]; 
    B = [OEI_Moment; OEI_Aileron];
    
    X = linsolve(A, B);
    forceRFInertOEI = X(1);
    forceRRInertOEI = X(2);
    
    SparReactInertLoadTot = total/9.81; 
    
    SparReactInertLoadTot(31) = SparReactInertLoadTot(31) - forceRFInertOEI;
    SparReactInertLoadTot(39) = SparReactInertLoadTot(39) - forceRRInertOEI;
    
 %   SparReactInertLoadTot(OEI_Position) = SparReactInertLoadTot(OEI_Position) + OEI_Thrust;
    
    SF_OEI = [];
    SF_OEI(1) = SparReactInertLoadTot(1);
    for i = 2:length(SparReactInertLoadTot)
        SF_OEI(i) = SF_OEI(i-1) + SparReactInertLoadTot(i);
    end

    figure
    plot(xDiscr, SF_OEI)

end