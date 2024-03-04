function tauFunc = tauFuncSim(alpha, beta, numSim)

% produces a simulated estimate for tau + 1/beta * log(lambda_tau + alpha)
% where tau is the natural closure time. 
% leverages the conditional uniformity of the cluster

    rho = alpha/beta;
    CDF = exp(-alpha/beta);

    tauFunc = 0;

    for n = 1 : numSim

        [N, CDF] = randBorel(rho, CDF);

        if N > 1

            park = randPF3(N-1);

            U = rand(1,N-1);
            Lambda = [0 sort(park - U)];

            tau = 0;
            for i = 2 : N
                tau = tau - log((i-1 - Lambda(i))/(i-1 - Lambda(i-1))) / beta;
            end

            tauFunc = tauFunc + tau + 1/beta*log(alpha*(N - Lambda(N)));

        else

            tauFunc = tauFunc + 1/beta*log(alpha);

        end

    end

    tauFunc = tauFunc / numSim;

