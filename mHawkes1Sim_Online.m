function [Qs, Ps, Pchecks] = mHawkes1Sim_Online(alpha, beta, lambda, polInd, polCo, polPow,  T, numSim)

%% code from CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% % By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov

    Ps = zeros(numSim, 1);
    Pchecks = zeros(numSim, 1);
    Qs = zeros(numSim, 1);
    

    for i = 1 : numSim

        [premCount, Q, N, premCheck] = mHawkes1Rep_Online(alpha, beta, lambda, polInd, polCo, polPow, T);

        Ps(i) = premCount/(N - Q);
        Qs(i) = Q;
        Pchecks(i) = premCheck/(N - Q);

    end

end
