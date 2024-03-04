function [ET, VT, P] = uhpAbsorbCheckSim(alpha, beta, polInd, checkVal, numSim)

    ET = 0;
    VT = 0;
    P = 0;

    for i = 1 : numSim

        if polInd == 1
            [time, prem] = uhpAbsorbCheckRep(alpha, beta, checkVal, Inf);
        elseif polInd == 2
            [time, prem] = uhpAbsorbCheckRep(alpha, beta, alpha, checkVal);
        end
        ET = ET + time/numSim;
        VT = VT + time^2/numSim;
        P = P + prem;

    end
    
    VT = VT - ET^2;
    P = P / numSim;

end