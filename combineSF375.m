function [total,comboSF] = combineSF375(xDiscr,InertialLoads,aero)

   
    totalLoad = sum(InertialLoads)+sum(aero);

    for i = 1:length(xDiscr)
        inertialMoments(i) = (InertialLoads(i)+aero(i)) * xDiscr(i);
    end

    totalMoment = sum(inertialMoments);

%    fuselageCG = totalMoment/totalLoad;

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
        total(i) = InertialLoads(i)*3.75 + aero(i) + sparReaction(i);
    end
    

    comboSF = total(1);
    for i = (2:length(total))
        comboSF(i) = comboSF(i-1) + total(i);
    end


end