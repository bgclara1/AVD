function [total,comboSF] = combineSFneg15(xDiscr,InertialLoads,LHT15) % effectively the n=3.75 case
   
    Load = InertialLoads*-1.5;
    Load(71) = Load(71) + LHT15/2*-1;
    Load(75) = Load(75) + LHT15/2*-1;
    totalLoad = sum(Load);

    for i = 1:length(xDiscr)
        inertialMoments(i) = Load(i) * xDiscr(i);
    end

    totalMoment = sum(inertialMoments);

    A = [ 31	39 ;
            1	1];
    B = [totalMoment; totalLoad];
    
    X = linsolve(A,B);
    forceRF = X(1);
    forceRR = X(2);


    sparReaction = zeros(1,80);
    sparReaction(31) = forceRF*-1;
    sparReaction(39) = forceRR*-1;

    for i = 1:length(xDiscr)
        total(i) = Load(i)  + sparReaction(i);
    end
    
    

    comboSF = total(1);
    for i = (2:length(total))
        comboSF(i) = comboSF(i-1) + total(i);
    end


end