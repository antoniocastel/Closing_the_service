function [N, CDFlocal] = randBorel(rho, CDFlocal)

%% code from CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% % By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov


%% Originally from 2024 “Conditional uniformity and Hawkes processes”. Mathematics of Operations Research 49(1):40–57 https: //doi.org/10.1287/moor.2022.1348.
    %% by Andrew Daw


    % method 1 - direct from CDF:

%     U0 = rand();
%     pN = exp(-rho);
%     N = 1;
%     while U0 > pN
%         N = N + 1;
%         pN = pN + exp(-rho*N - gammaln(N+1) + (N-1)*log(rho*N));
%     end

    % method 2 - Poisson random walk:
    
%     X = 1;
%     N = 1;
%     while X > 0
%         newX = poissrnd(rho);
%         N = N + newX;
%         X = X + newX - 1;
%     end

    % method 3 - CDF with adaptive precompute:

    U0 = rand();
    if CDFlocal(end) > U0
        
        N = find(U0 <= CDFlocal, 1);
        
    else
        
        N = numel(CDFlocal);
        while U0 > CDFlocal(N)
            N = N + 1;
            CDFlocal(N) = CDFlocal(N-1) + exp(-rho*N - gammaln(N+1) + (N-1)*log(rho*N));
        end
        
    end
    
end