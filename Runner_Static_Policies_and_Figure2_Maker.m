
%% Code For Static Polcies  
% and to generate Figure 1 of CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov



alphaOverBeta=[0.8,0.9,0.95,0.99];

betaVec      =  [4;8;16];

numBetas     = numel(betaVec);
alphaMat     = betaVec.*alphaOverBeta;
numAlphas    = size(alphaMat,2);
parCount     = numel(alphaMat);

lambda = 1;

numSimCheck = 2^20;

n=100;
resEQ1     = NaN(n,parCount);
resES1     = NaN(n,parCount);
resP1      = NaN(n,parCount);

resEQ2     = NaN(n,parCount);
resES2     = NaN(n,parCount);
resP2      = NaN(n,parCount);

colum_no=1;

for i = 1 : numBetas
    beta = betaVec(i);

    for j = 1: numAlphas
        alpha = alphaMat(i,j);
        ratio = alphaOverBeta(j);


        % finding max delta via simple binary search 
        tic
        tauFunc = tauFuncSim(alpha, beta, 2^22);
        toc
        maxDeltaL = max(0, alpha - exp(-beta*(1 - tauFunc)));
        maxDeltaU = alpha*(1- exp(-beta));
        nearStab = false;
        minP = 0.9999;
        while ~nearStab

            maxDelta = (maxDeltaU + maxDeltaL)/2;
            tic
            ET = uhpAbsorbCheckSim(alpha, beta, 1, maxDelta, numSimCheck);
            toc
            if ET < 1 && ET > minP
                nearStab = true;
            elseif ET > minP
                maxDeltaU = maxDelta;
            else
                maxDeltaL = maxDelta;
            end

        end

        ptStep = 0.01*(1 - exp(-alpha/beta));
        maxPt = (ptStep^2) * floor((1 - exp(-alpha/beta))/(ptStep^2));
        minPt = ptStep * ceil((1 - exp(-(alpha-maxDelta)/beta))/ptStep);
        ptVec  = minPt : ptStep : maxPt;
        num_pt  = numel(ptVec);




        %% static policy delta policy
        polInd = 1;
        for jj=1:num_pt

            de = alpha + beta*log(1-ptVec(jj));

            tic
            [ET, VT, Pc] = uhpAbsorbCheckSim(alpha, beta, polInd, de, numSimCheck);
            toc
            eqc = (ET^2 + VT)/(2*(1 - ET));
            resEQ1(jj,colum_no) = eqc;
            resP1(jj,colum_no)  =  Pc;
            resES1(jj,colum_no) = ET;


        end


        %% static epsilon policy
        polInd = 2;
        epsStep = 0.025;
        eps = 0.05;
        esn = 0;

        stabBound = max(resES1(:,colum_no));
      
        meanTime = 0;
        counter = 1;
        while meanTime < stabBound
            tic
            [ET, VT, Pc] = uhpAbsorbCheckSim(alpha, beta, polInd, eps, numSimCheck);
            toc
            eqc = (ET^2 + VT)/(2*(1 - ET));


            resEQ2(counter,colum_no) = eqc;
            resP2(counter,colum_no)  = Pc;
            resES2(counter,colum_no) = ET;

            meanTime = ET;
            counter=counter+1;
            eps = eps + epsStep;
        end

        colum_no=        colum_no+1;
    end
end

%% Save results 

save('StaticPolicies_res')



%% Plot


colorsArrayBlue = {
    [0.0314, 0.2706, 0.5804],
    [0.2588, 0.5725, 0.7765], 
    [0.4196, 0.6824, 0.8392], 
    [0.6196, 0.7922, 0.8824], 
    [0.7765, 0.8588, 0.9373], 
    [0.9373, 0.9529, 1.0000]
};




colorsArrayYellow = {
    [0.6000, 0.2039, 0.0157], 
    [0.9961, 0.6000, 0.1608],
    [0.9961, 0.7686, 0.3098], 
    [0.9961, 0.8902, 0.5686], 
    [1.0000, 0.9686, 0.7373], 
    [1.0000, 1.0000, 0.8980]
};

colorsArray = {colorsArrayBlue,  
    colorsArrayYellow
};


policyMarkers = {'o', '*', 'square', '^','|','diamond'}; % Different markers for each policy



for j = 1: numAlphas

    ratio = alphaOverBeta(j);
    colum_no=        j;

    %For Figures
    legendEntries = [];
    legendLabels = {};
    ratioLegendEntries = [];
    ratioLegendLabels = {};

    figure(j); % Create a new figure for each ratio
    hold on;


    for i= 1 : numBetas
        alpha = alphaMat(i,j);

        beta = betaVec(i);

        polInd = 1;


        % Select the color array for the current policy index
        currentColorArray = colorsArray{polInd};

        % Compute color index within the array
        colorIndex = min(i, length(currentColorArray)); 
        color = currentColorArray{colorIndex};
         


        name= 'Level';

            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', [name, ' (\beta= ' ,  num2str(beta), ')'], 'MarkerEdgeColor', color,'MarkerFaceColor',color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = [name, ' (\beta= ' ,  num2str(beta), ')'];


        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)),  'Color', color,'LineWidth',1.5);


        polInd = 2;
         % Select the color array for the current policy index
        currentColorArray = colorsArray{polInd};

        % Compute color index within the array, adjust logic as needed
        colorIndex = min(i, length(currentColorArray)); 
        color = currentColorArray{colorIndex};

        name= 'Time';
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', [name, ' (\beta= ' ,  num2str(beta), ')'], 'MarkerEdgeColor', color,'MarkerFaceColor',color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = [name, ' (\beta= ' ,  num2str(beta), ')'];
        plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.5);

        
        colum_no=        colum_no+4;

    end


  


    hold off;
    legend(legendEntries,  legendLabels,'Location', 'southwest');
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Probability (Active Customers)' )
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 16);
    xlim([10^-3.2, 10^3]); % Set x-axis limits

    ylim([10^-3.05, 10^0]); % Set y-axis limits

    name =  sprintf('V8_Static_ParetoF_rat%.2f.fig', ratio);
    name2 = sprintf('V8_Static_ParetoF_rat%.2f.jpg',  ratio);
    savefig(name);
    saveas(gcf, name2);
    name2 = sprintf('V8_Static_ParetoF_rat%.2f.eps',  ratio);
    saveas(gcf, name2);
    close;


end
