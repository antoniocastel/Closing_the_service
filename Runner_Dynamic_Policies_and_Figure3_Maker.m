

%% Code For Dynamic Polcies  % and to generate Figure 3 of CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov


makePlots = true;

%alphaOverBeta
rhoVec = [0.8,0.9,0.95,0.99];
betaVec = [4; 8; 16];

colorMat = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

numRho = numel(rhoVec);
numBeta = numel(betaVec);

numSim = 2^14; %14
T = 2^13; %12

lambda = 1;

polPowVec = [0.5; 1; 2];
numPolPow = numel(polPowVec);
numPolCo = 20;
polCoMat(1, 1:numPolCo) = 2.^(-(0:1:(numPolCo-1)));
polCoMat(2, 1:numPolCo) = (3/2).^((-1):1:(numPolCo-2));

markerMat = ['<', '^', '>'; '+', '*', 'x'];

for rhoInd = 1 : numRho

    rho = rhoVec(rhoInd);

    if makePlots

        figure
        hold on
        set(gca, 'yscale', 'log')
        set(gca, 'xscale', 'log')

    end

    for betaInd = 1 : numBeta

        beta = betaVec(betaInd);
        alpha = rho*beta;

        for polPowInd = 1 : numPolPow

            polPow = polPowVec(polPowInd);

            for polCoInd = 1 : numPolCo
    
                tic
    
                for polInd = 1 : 2 % level and then time
    
                    polCo = polCoMat(polInd, polCoInd);
                    if polInd == 1
                        polCo = beta*polCo^polPow;
                    else
                        polCo = 1/beta*polCo^polPow;
                    end
    
                    [Qs, Ps, ~] = mHawkes1Sim_Online(alpha, beta, lambda, polInd, polCo, polPow, T, numSim);
        
                    if polInd == 1
        
                        EQ1(polCoInd, polPowInd, betaInd, rhoInd) = mean( max(Qs-1,0) );
                        P1(polCoInd, polPowInd, betaInd, rhoInd) = mean( Ps );
        
                    else
        
                        EQ2(polCoInd, polPowInd, betaInd, rhoInd) = mean( max(Qs-1,0) );
                        P2(polCoInd, polPowInd, betaInd, rhoInd) = mean( Ps );
        
                    end
    
                end
    
                clockCheck = [toc, rhoInd, betaInd, polPowInd, polCoInd]
    
            end


            if makePlots

                hold on

                plotQ1s = EQ1(:,polPowInd,betaInd,rhoInd);
                plotP1s = P1(:,polPowInd,betaInd,rhoInd);
                plotQ2s = EQ2(:,polPowInd,betaInd,rhoInd);
                plotP2s = P2(:,polPowInd,betaInd,rhoInd);
    
                plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['Level (\beta = ', num2str(beta), ', k = ', num2str(polPow), ')'], 'MarkerFaceColor', colorMat(betaInd,:))
                
                plot(plotQ2s(plotQ2s < T/8), plotP2s(plotQ2s < T/8), 'linestyle', '-.', 'marker', markerMat(2, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['Time (\beta = ', num2str(beta), ', k = ', num2str(polPow), ')'])
    
            end

        end

    end

    if makePlots

        legend('location', 'southwest')

    end
    
end

save('DynamicPolicies_res')


%% Plot
% load Static Results 

load('StaticPolicies_res')



colorMat = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
    markerMat = ['<', '^', '>'; '+', '*', 'x'];
    Indexes=[[1,5,9];[2,6,10];[3,7,11];[4,8,12]];

for rhoInd = 1 : numRho

    rho = rhoVec(rhoInd);

    indexes = Indexes(rhoInd,:);



    figure
    hold on
    set(gca, 'yscale', 'log')
    set(gca, 'xscale', 'log')
    hold on

    column = 1;

    legendHandles = [];
    betaInd = 1;
    h =         plot(resEQ1(:, indexes(betaInd)), resP1(:, indexes(betaInd))./ (1 - exp(-alpha /beta)),...
        'Color', [0.2,0.2,0.2], 'DisplayName', 'Static' , 'LineWidth', 1);
    legendHandles(end+1) = h;

    for betaInd = 2 : numBeta

        beta = betaVec(betaInd);
        alpha = rho*beta;
        
        plot(resEQ1(:, indexes(betaInd)), resP1(:, indexes(betaInd))./ (1 - exp(-alpha /beta)),...
            'Color', [0.2,0.2,0.2] , 'LineWidth', 1);

    end

    betaInd = 1 ;
    beta = betaVec(betaInd);
    alpha = rho*beta;



    polPowInd = 1;
    polPow = polPowVec(polPowInd);

    plotQ1s = EQ1(:,polPowInd,betaInd,rhoInd);
    plotP1s = P1(:,polPowInd,betaInd,rhoInd)./(1 - exp(-alpha /beta));
    plotQ2s = EQ2(:,polPowInd,betaInd,rhoInd);
    plotP2s = P2(:,polPowInd,betaInd,rhoInd)./(1 - exp(-alpha /beta));


    h = plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\ell(q)=\theta q^2$, ($\beta=$', num2str(beta),')'],  'MarkerFaceColor', colorMat(betaInd,:));
    legendHandles(end+1) = h;



    h2 = plot(plotQ2s(plotQ2s < T/8), plotP2s(plotQ2s < T/8), 'linestyle', '-.', 'marker', markerMat(2, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\epsilon (q)=\Delta/q$, ($\beta=$', num2str(beta),')']);
    legendHandles(end+1) = h2;

    polPowInd = 2;
    h = plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\ell(q)=\theta q$, $\beta=$', num2str(beta)],  'MarkerFaceColor', colorMat(betaInd,:));
    legendHandles(end+1) = h;
    h2 = plot(plotQ2s(plotQ2s < T/8), plotP2s(plotQ2s < T/8), 'linestyle', '-.', 'marker', markerMat(2, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\epsilon (q)=\Delta/q^{2}$, ($\beta=$', num2str(beta),')']);
    legendHandles(end+1) = h2; %
    polPowInd = 3;
    h = plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\ell(q)=\theta \sqrt{q}$, $\beta=$', num2str(beta)],  'MarkerFaceColor', colorMat(betaInd,:));
    legendHandles(end+1) = h;
    h2 = plot(plotQ2s(plotQ2s < T/8), plotP2s(plotQ2s < T/8), 'linestyle', '-.', 'marker', markerMat(2, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['$\epsilon (q)=\Delta/\sqrt{q}$, ($\beta=$', num2str(beta),')']);
    legendHandles(end+1) = h2; %


    for betaInd = 2 : numBeta
        beta = betaVec(betaInd);
        alpha = rho*beta;


        for polPowInd = 1 : numPolPow
            polPow = polPowVec(polPowInd);

            plotQ1s = EQ1(:,polPowInd,betaInd,rhoInd);
            plotP1s = P1(:,polPowInd,betaInd,rhoInd)./(1 - exp(-alpha /beta));
            plotQ2s = EQ2(:,polPowInd,betaInd,rhoInd);
            plotP2s = P2(:,polPowInd,betaInd,rhoInd)./(1 - exp(-alpha /beta));


            if polPowInd==1
                h = plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:), 'DisplayName', ['($\beta=$', num2str(beta),')'],  'MarkerFaceColor', colorMat(betaInd,:));
                legendHandles(end+1) = h; % Add handle to legend array

            else
                plot(plotQ1s(plotQ1s < T/8), plotP1s(plotQ1s < T/8), 'marker', markerMat(1, polPowInd), 'Color', colorMat(betaInd,:),  'MarkerFaceColor', colorMat(betaInd,:))
            end


            plot(plotQ2s(plotQ2s < T/8), plotP2s(plotQ2s < T/8), 'linestyle', '-.', 'marker', markerMat(2, polPowInd), 'Color', colorMat(betaInd,:))

        end

    end

    legend(legendHandles, 'Location', 'SouthWest', 'Interpreter', 'latex', FontWeight='bold');

    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Probability (Active Customers)' )
    xlim([10^-3.2, 10^3]); % Set x-axis limits

    ylim([10^-7, 10^0]); % Set y-axis limits
    set(gca, 'FontSize', 16);

    text(10^-3, 0.1, '$\beta=16$', 'FontSize', 16, 'VerticalAlignment', 'baseline', Interpreter='latex', FontWeight='bold');
    text(10^-1, 0.2, '$\beta=8$', 'FontSize', 16, 'VerticalAlignment', 'baseline', Interpreter='latex', FontWeight='bold');
    text(10^1, 0.1, '$\beta=4$', 'FontSize', 16, 'VerticalAlignment', 'baseline', Interpreter='latex', FontWeight='bold');

    hold off
    name =  sprintf('V4LongRun_Online_ParetoF_rat%.2f.fig', rho);
    name2 = sprintf('V4LongRun_Online_ParetoF_rat%.2f.jpg',  rho);
    savefig(name);
    saveas(gcf, name2);
    name2 = sprintf('V4LongRun_Online_ParetoF_rat%.2f.eps',  rho);
    saveas(gcf, name2);
    close;

end

