function [premCount, Q, N, premCheck] = mHawkes1Rep_Online(alpha, beta, lambda, polInd, polCo, polPow, T)

%% code from CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% % By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov

    Q = 0;
    time = 0;
    mu = 0;
    N = 0;
    premCount = 0;
    premCheck = 0;
    last = 0;

    while time < T

        % time to next arrival
        S1 = - log(rand())/lambda;

        % time until next message and time until closure
        if Q > 0
            % time to next message (if happens) via Dassios-Zhao
            % check if happens (D > 0)
            D = 1 + beta*log(rand())/mu;
            if D < 0
                % time to next message is infinite
                S2 = Inf;
            else
                % time to next message
                S2 = -log(D)/beta;
            end

            % time until systematic closure, depending on policy
            if polInd == 1 
                % online level policy 
                level = polCo * (Q-1)^polPow;
                S3 = max(0, log(mu/level)/beta);
            elseif polInd == 2
                % online time policy 
                eps = polCo / (Q-1)^polPow;
                S3 = max(last + eps - time, 0);
            end

        else
            % no ongoing service
            S2 = Inf;
            S3 = Inf;
        end

        % time until end of sim
        S4 = T - time;

        if S1 < min([S2 S3 S4])
            % arrival next
            N = N + 1;
            Q = Q + 1;
            time = time + S1;
            if Q > 1
                mu = mu*exp(-beta*S1);
            else
                mu = alpha;
                last = time;
            end

        elseif S2 < min([S1 S3 S4])
            % message next
            time = time + S2;
            mu = mu*exp(-beta*S2) + alpha;
            last = time;

        elseif S3 < min([S1 S2 S4])
            % end of service
            time = time + S3;
            if S2 < Inf 
                premCount = premCount + 1;
            end


            %%%
            Q = Q - 1;
            if Q > 0
                mu = alpha;
                last = time;
            else
                mu = 0;
            end

        else
            % end of simulation
            time = T;
            mu = mu*exp(-beta*S4);
            
        end

    end

end