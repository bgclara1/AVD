function [total,comboSF] = combineSF375(xDiscr,InertialLoads,LHT)

    Load = InertialLoads*3.75;
    Load(71) = Load(71)+ LHT/2;
    Load(75) = Load(75)+ LHT/2;
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
    disp('3.75')
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