function [time, prem] = uhpAbsorbCheckRep(alpha, beta, delta, eps)


    time = 0;
    S = 0;
    prem = 0;
    S1 = 0;
    S2 = Inf;

    mu = alpha;
    
    while S1 < S2
       
        % S1 is next message time (if finite, per D)
        D = 1 + beta*log(rand)/mu;
        if D < 0
            S1 = Inf;
        else
            S1 = -log(D)/beta;
        end
        
        % S2 is closure time:
        S2 = min([log(mu/(alpha - delta))/beta, eps]);
        
        % S is next event time:
        S = min([S1, S2]);
        
        time = time + S;
        mu = mu*exp(-beta*S);
        
        if S1 < S2
            mu = mu + alpha;
        elseif D > 0
            prem = 1;
        end
        
    end
    
end
    