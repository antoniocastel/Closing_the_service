function [Qs, Ps, ES] = online_mHawkes1Sim(alpha, beta, lambda, polInd, de, T, numSim,k)

    Ps = zeros(numSim, 1);
    Qs = zeros(numSim, 1);
    ES = zeros(numSim, 1);
    ESdenom =zeros(numSim, 1);
    %ES = 0;
    %ESdenom = 0;
        if polInd == 4
            %tic
            %index 1 is 0 in queue
            falseQ = 1;
            Lambert_0_Vec = [];
            while de*(falseQ-1) <= beta*exp(-1)
                Lambert_0_Vec = [Lambert_0_Vec; lambertw(0, -de * (falseQ-1) / beta)];
                if ~isreal(Lambert_0_Vec(falseQ))
                    out = 'error -- complex value from lambert w'
                    return;
                end
                falseQ = falseQ + 1;
            end
            %toc 
        else
            Lambert_0_Vec=0;
        end 
    
       
    parfor i = 1 : numSim
        % if  polInd == 4
        %     [premCount, Q, N, totalS] = online_mHawkes1Rep_with_tdelta(alpha, beta, lambda, polInd, de, T, Lambert_0_Vec);
        % else 
        %     [premCount, Q, N, totalS] = online_mHawkes1Rep_with_tdelta(alpha, beta, lambda, polInd, de, T);
        % end 

        [premCount, Q, N, totalS] = online_mHawkes1Rep_with_tdelta(alpha, beta, lambda, polInd, de, T, Lambert_0_Vec,k );

        Ps(i)      = premCount/(N - Q);
        Qs(i)      = Q;
        ES(i)      = totalS;
        ESdenom(i) = N - Q;
        % ES = ES + totalS;
        % ESdenom = ESdenom + N - Q;

    end
    
    ES      =sum(ES);
    ESdenom =sum(ESdenom);
    %Pcalc = 1 - exp(-(alpha - de)/beta);

    if ES~=0
        ES = ES / ESdenom;
    end    

end