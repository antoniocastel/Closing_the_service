addpath('/Users/antoniomoshe/Library/CloudStorage/Dropbox/Chat systems - Service time modeling/ClosingTheService/Simulations/M-Hawkes-1')


alphaOverBeta=[0.5,0.75,0.9];

betaVec      =  [2;4;6;8];
numBetas     = numel(betaVec);
alphaMat     = betaVec.*alphaOverBeta;
numAlphas    = size(alphaMat,2);
parCount     = numel(alphaMat);

lambda = 1;
%sensitivity parameter to waiting customers
% Theta     = 0.1:0.1:1;
% numThetas = numel(Theta); 
minTheta = 0.001;
minC     = 0.001;
minK     = 0.001;



numSimCheck = 2^20;
numSim      = 2^14;
T           = 2^10;

%polInd = 3;
%polInd = 4;

%vectors to save results
n=100;
resEQ1     = NaN(n,parCount);
resES1     = NaN(n,parCount);
resP1      = NaN(n,parCount);

resEQ2     = NaN(n,parCount);
resES2     = NaN(n,parCount);
resP2      = NaN(n,parCount);
n = ceil(log2(max(alphaMat,[],'all')/minTheta));
resEQ3     = NaN(n,parCount);
resES3     = NaN(n,parCount);
resP3      = NaN(n,parCount);

resEQ4     = NaN(n,parCount);
resES4     = NaN(n,parCount);
resP4      = NaN(n,parCount);

resEQ5     = NaN(n,parCount);
resES5     = NaN(n,parCount);
resP5      = NaN(n,parCount);

resEQ6     = NaN(n,parCount);
resES6     = NaN(n,parCount);
resP6      = NaN(n,parCount);


minK     = 1;

n = ceil(log2(max(alphaMat,[],'all')/minK));
n2 = ceil(log2(max(alphaMat,[],'all')/minK));

resEQ7     = NaN(n*n2,parCount);
resES7     = NaN(n*n2,parCount);
resP7      = NaN(n*n2,parCount);


colum_no=1;


%   policyMarkers = {'o', 'x', 'square', '^','.'}; % Different markers for each policy
%policyMarkers = {'o', '*', 'square', '^','|'}; % Different markers for each policy
policyMarkers = {'o', '*', 'square', '^','|','diamond'}; % Different markers for each policy

baseRatioColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880],[1, 1, 0]}; % Blue, Orange, Green in RGB



for i = 1 : numBetas
    beta = betaVec(i);
   
    %For Figures
    legendEntries = []; 
    legendLabels = {}; 
    ratioLegendEntries = []; 
    ratioLegendLabels = {}; 

    figure(i); % Create a new figure for each beta value
    hold on; 
    ylim([0, 1]);

    for j = 1: numAlphas
        alpha = alphaMat(i,j);
        ratio = alphaOverBeta(j);


        % finding max delta via simple binary search (stop at first in (0.97, 1))
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

        baseColor = baseRatioColors{j};
        shadeFactor = linspace(0.8, 1.2, 5); % shade intensity
        
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static Level Pol';
        
        if j == 1
            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ1(:, colum_no),  resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        %% static epsilon policy
        polInd = 2;
        epsStep = 0.025;
        eps = 0.05;
        esn = 0;
        %if numel(ES) > 0
        %    stabBound = max(ES); % could also just do 1, but this is perhaps more fair
        %else
            stabBound = max(resES1(:,colum_no));
        %end
        meanTime = 0;
        counter = 1; 
        while meanTime < stabBound
            tic
            [ET, VT, Pc] = uhpAbsorbCheckSim(alpha, beta, polInd, eps, numSimCheck);
            toc
            eqc = (ET^2 + VT)/(2*(1 - ET));
            % EQ2check = [EQ2check; eqc];
            % P2check = [P2check; Pc];
            
            resEQ2(counter,colum_no) = eqc;
            resP2(counter,colum_no)  = Pc;
            resES2(counter,colum_no) = ET;

            meanTime = ET;
            counter=counter+1;
            eps = eps + epsStep;
        end
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static Time Pol';
        if j == 1
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ2(:, colum_no),  resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        h2= plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'DisplayName', ['\alpha/\beta= ' ,  num2str(ratio)], 'Color', color,'LineWidth',1.5);
        ratioLegendEntries = [ratioLegendEntries, h2];
        ratioLegendLabels{end+1} = ['\alpha/\beta= ' ,  num2str(ratio)];

        %% Online theta Policy
        polInd = 3; 
        n = ceil(log2(alpha/minTheta));
        Thetas = alpha * 2.^(-linspace(0, n, n+1));

        
        % numThetas = numel(Theta);
        %Set the sensitivity parameter theta
        for ii = 1: n
        %for ii = 1: numThetas
            th=Thetas(1,ii);

            tic
            [qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, th, T, numSim);
            toc

            resEQ3(ii,colum_no) = mean( max(qs-1,0) );
            resP3 (ii,colum_no) = mean( ps );
            resES3(ii,colum_no) = esn;
            
        end 
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta-Pol';
        if j == 1
            h = scatter(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ3(:, colum_no),  resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);

        
        %% Online c Policy
        polInd = 4; 
        n = ceil(log2(alpha/minC));
        Cs = alpha * 2.^(-linspace(0, n, n+1));
        for ii = 1: n
        %for ii = 1: numThetas
            c=Cs(1,ii);

            tic
            [qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, c, T, numSim);
            toc

            resEQ4(ii,colum_no) = mean( max(qs-1,0) );
            resP4 (ii,colum_no) = mean( ps );
            resES4(ii,colum_no) = esn;
            
        end 
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online c-Pol';
        if j == 1
            h = scatter(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ4(:, colum_no),  resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        %% Online theta q^2 Policy
        polInd =5; 
        n = ceil(log2(alpha/minTheta));
        Thetas = alpha * 2.^(-linspace(0, n, n+1));
        for ii = 1: n
        %for ii = 1: numThetas
            th=Thetas(1,ii);

            tic
            [qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, th, T, numSim);
            toc

            resEQ5(ii,colum_no) = mean( max(qs-1,0) );
            resP5 (ii,colum_no) = mean( ps );
            resES5(ii,colum_no) = esn;
            
        end 
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta*q^2';
        if j == 1
            h = scatter(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ5(:, colum_no),  resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);
        
        %% Online theta*q^(1/2)

        polInd =6; 
        n = ceil(log2(alpha/minTheta));
        Thetas = alpha * 2.^(-linspace(0, n, n+1));
        for ii = 1: n
        %for ii = 1: numThetas
            th=Thetas(1,ii);

            tic
            [qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, th, T, numSim);
            toc

            resEQ6(ii,colum_no) = mean( max(qs-1,0) );
            resP6 (ii,colum_no) = mean( ps );
            resES6(ii,colum_no) = esn;
            
        end 

         
        colum_no=        colum_no+1;

         

    end 
    hold off;
    legend([legendEntries, ratioLegendEntries], [legendLabels, ratioLegendLabels],'Location', 'best'); % Combine policy and ratio legends
    title(['Pareto Frontiers for \beta = ', num2str(beta)]);
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Rate for Active Customers')
    
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 14);
    name =  sprintf('Stc_Onl_ParetoF_bet%.1f_lmb%.1f.fig', beta, lambda);
    name2 = sprintf('Stc_Onl_ParetoF_bet%.1f_lmb%.1f.jpg',  beta, lambda);
    savefig(name);
    saveas(gcf, name2);
    close;
end


%% 
save('res_5_Policies.mat')
save('res_6_Policies.mat')

%% Policy 7 (online) Virtual Abandonment is distributed Gamma 


Thetas=1:1:6;
Ks    =2:3;
n=length(Ks)*length(Thetas);
resEQ7     = NaN(n,parCount);
resES7     = NaN(n,parCount);
resP7      = NaN(n,parCount);



colum_no=1;

for i =4:numBetas
%for i = 1 : numBetas
    beta = betaVec(i);

    for j = 1: numAlphas
        ratio = alphaOverBeta(j);
        alpha = alphaMat(i,j);
        polInd =7; 
        
        %n = ceil(log2((alpha)/minK));
        %Thetas = alpha * 2.^(-linspace(0, n, n+1));
        
        %n2 = ceil(log2((alpha)/minK));
        %Ks = alpha * 2.^(-linspace(0, n, n+1));
        cnt=1;
        %for ii = 4: length(Thetas)
        for ii = 1: length(Thetas)
            th=Thetas(1,ii);
           for ll = 1: length(Ks)
               k=Ks(1,ll);

               tic
               %[qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, th, T, numSim,k);
               [qs, ps, esn] = online_mHawkes1Sim(alpha, beta, lambda, polInd, th, T, numSim,k);
               toc
              
               resEQ7(cnt,colum_no) = mean( max(qs-1,0) );
               resP7 (cnt,colum_no) = mean( ps );
               resES7(cnt,colum_no) = esn;

                cnt=cnt+1;
                save('res_7_i.mat','resEQ7','resP7','resES7')

           end
            
        end 

         
        colum_no=        colum_no+1;
    end
end


%% Save
save('res_7_Policies.mat')

% Sort each column 
sortedP7 = zeros(size(resP7));
sortedQ7 = zeros(size(resEQ7));

% Loop through each row
for i = 1:size(resP7,2)
    [sortedP7(:, i), sortIndices] = sort(resP7(:, i), 'descend');
    
    sortedEQ7(:, i) = resEQ7(sortIndices,i);
end


%% V3. coloring by Policy Type and Shading by Beta. 
%  Plot maker one Plot for each alpha/beta


baseRatioColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880],[0.4660, 0.6740, 0.1880],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840],[0.9290 0.6940 0.1250]}; % Blue, Orange, Green in RGB

policyMarkers = {'o', '*', 'square', '^','^','|','diamond'}; % Different markers for each policy

for j = 1: numAlphas

    ratio = alphaOverBeta(j);
    colum_no=        j;

    %For Figures
    legendEntries = [];
    legendLabels = {};
    ratioLegendEntries = [];
    ratioLegendLabels = {};

    figure(j); % Create a new figure for each beta value
    hold on;
    ylim([0, 1]);
    shadeFactor = linspace(0.8, 1.2, 4); % shade intensity by beta

    for i = 1 : numBetas
        alpha = alphaMat(i,j);

        beta = betaVec(i);

        polInd = 1;

        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static level pol';

        if i == 1
            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ1(:, colum_no),  resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        polInd = 2;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static time pol.';
        if i == 1
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ2(:, colum_no),  resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end

        h2= plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'DisplayName', ['\beta= ' ,  num2str(beta)], 'Color', color,'LineWidth',1.2);
        ratioLegendEntries = [ratioLegendEntries, h2];
        ratioLegendLabels{end+1} =  ['\beta= ' ,  num2str(beta)];

        polInd = 3;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta-pol.';
        if i == 1
            h = scatter(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ3(:, colum_no),  resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);


        % polInd = 4;
        % baseColor = baseRatioColors{polInd};
        % 
        % color = baseColor * shadeFactor(i);
        % color = min(max(color, 0), 1); % Ensure color values are within valid range
        % name= 'Online c-pol.';
        % if i == 1
        %     h = scatter(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
        %     legendEntries = [legendEntries, h];
        %     legendLabels{end+1} = name;
        % else
        %     scatter(resEQ4(:, colum_no),  resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        % end
        % plot(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);

        polInd =5;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^2-pol.';
        if i == 1
            h = scatter(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ5(:, colum_no),  resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);
        


        polInd =6;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^{0.5}-pol.';
        if i == 1
            h = scatter(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ6(:, colum_no),  resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);

        polInd =7;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \Gamma-pol.';
        if i == 1
            h = scatter(sortedEQ7(:, colum_no), sortedP7(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(sortedEQ7(:, colum_no),  sortedP7(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(sortedEQ7(:, colum_no), sortedP7(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);



        colum_no=        colum_no+3;

    end
    hold off;
    legend([legendEntries, ratioLegendEntries], [legendLabels, ratioLegendLabels],'Location', 'best'); % Combine policy and ratio legends
    title(['Pareto Frontiers for \alpha/\beta = ', num2str(ratio)]);
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Rate for Active Customers')
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 14);
    name =  sprintf('V2_Stc_Onl7_ParetoF_rat%.1f_lmb%.1f.fig', ratio, lambda);
    name2 = sprintf('V2_Stc_Onl7_ParetoF_rat%.1f_lmb%.1f.jpg',  ratio, lambda);
    savefig(name);
    saveas(gcf, name2);
    close;
end










%% Plot maker one plot for each Beta

baseRatioColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880],[1, 1, 0]}; % Blue, Orange, Green in RGB

policyMarkers = {'o', '*', 'square', '^','|'}; % Different markers for each policy
colum_no=1;
for i = 1 : numBetas
    beta = betaVec(i);
    %For Figures
    legendEntries = []; 
    legendLabels = {}; 
    ratioLegendEntries = []; 
    ratioLegendLabels = {}; 

    figure(i); % Create a new figure for each beta value
    hold on; 
    ylim([0, 1]);

    for j = 1: numAlphas
        alpha = alphaMat(i,j);
        ratio = alphaOverBeta(j);

%%
        polInd = 1;
        baseColor = baseRatioColors{j};
        shadeFactor = linspace(0.8, 1.2, 5); % shade intensity

        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static level pol.';
        if j == 1
            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ1(:, colum_no),  resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);
%%
        polInd = 2;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static time pol.';
        if j == 1
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ2(:, colum_no),  resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        h2= plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'DisplayName', ['\alpha/\beta= ' ,  num2str(ratio)], 'Color', color,'LineWidth',1.5);
        ratioLegendEntries = [ratioLegendEntries, h2];
        ratioLegendLabels{end+1} = ['\alpha/\beta= ' ,  num2str(ratio)];

%%
        polInd = 3;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta-pol.';
        if j == 1
            h = scatter(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ3(:, colum_no),  resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);
%%
        polInd = 4;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online c-pol.';
        if j == 1
            h = scatter(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ4(:, colum_no),  resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);

        polInd =5;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^2-pol.';
        if j == 1
            h = scatter(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ5(:, colum_no),  resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);
        colum_no=        colum_no+1;



    end
    hold off;
    legend([legendEntries, ratioLegendEntries], [legendLabels, ratioLegendLabels],'Location', 'best'); % Combine policy and ratio legends
    title(['Pareto Frontiers for \beta = ', num2str(beta)]);
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Rate for Active Customers')
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 14);
    name =  sprintf('B_Stc_Onl_ParetoF_bet%.1f_lmb%.1f.fig', beta, lambda);
    name2 = sprintf('B_Stc_Onl_ParetoF_bet%.1f_lmb%.1f.jpg',  beta, lambda);
    savefig(name);
    saveas(gcf, name2);
    close;
end



%% Plot maker one Plot for each alpha/beta


baseRatioColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880],[0.9290 0.6940 0.1250]}; % Blue, Orange, Green in RGB

policyMarkers = {'o', '*', 'square', '^','|','diamond'}; % Different markers for each policy

for j = 1: numAlphas

    ratio = alphaOverBeta(j);
    colum_no=        j;

    %For Figures
    legendEntries = [];
    legendLabels = {};
    ratioLegendEntries = [];
    ratioLegendLabels = {};

    figure(j); % Create a new figure for each beta value
    hold on;
    ylim([0, 1]);

    for i = 1 : numBetas
        alpha = alphaMat(i,j);

        beta = betaVec(i);

        polInd = 1;

        baseColor = baseRatioColors{i};
        shadeFactor = linspace(0.8, 1.2, 6); % shade intensity

        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static level pol';

        if i == 1
            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ1(:, colum_no),  resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        polInd = 2;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static time pol.';
        if i == 1
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ2(:, colum_no),  resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        h2= plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'DisplayName', ['\beta= ' ,  num2str(beta)], 'Color', color,'LineWidth',1.5);
        ratioLegendEntries = [ratioLegendEntries, h2];
        ratioLegendLabels{end+1} =  ['\beta= ' ,  num2str(beta)];

        polInd = 3;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta-pol.';
        if i == 1
            h = scatter(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ3(:, colum_no),  resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        polInd = 4;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online c-pol.';
        if i == 1
            h = scatter(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ4(:, colum_no),  resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);

        polInd =5;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^2-pol.';
        if i == 1
            h = scatter(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ5(:, colum_no),  resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);
        


        polInd =6;
        color = baseColor * shadeFactor(polInd);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^{0.5}-pol.';
        if i == 1
            h = scatter(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ6(:, colum_no),  resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end
        plot(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        colum_no=        colum_no+3;

    end
    hold off;
    legend([legendEntries, ratioLegendEntries], [legendLabels, ratioLegendLabels],'Location', 'best'); % Combine policy and ratio legends
    title(['Pareto Frontiers for \alpha/\beta = ', num2str(ratio)]);
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Rate for Active Customers')
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 14);
    name =  sprintf('Stc_Onl6_ParetoF_rat%.1f_lmb%.1f.fig', ratio, lambda);
    name2 = sprintf('Stc_Onl6_ParetoF_rat%.1f_lmb%.1f.jpg',  ratio, lambda);
    savefig(name);
    saveas(gcf, name2);
    close;
end



%% V2. coloring by Policy Type and Shading by Beta. 
%  Plot maker one Plot for each alpha/beta


baseRatioColors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840]}; % Blue, Orange, Green in RGB

%policyMarkers = {'o', '*', 'square', '^','|','diamond'}; % Different markers for each policy

for j = 1: numAlphas

    ratio = alphaOverBeta(j);
    colum_no=        j;

    %For Figures
    legendEntries = [];
    legendLabels = {};
    ratioLegendEntries = [];
    ratioLegendLabels = {};

    figure(j); % Create a new figure for each beta value
    hold on;
    ylim([0, 1]);
    shadeFactor = linspace(0.8, 1.2, 4); % shade intensity by beta

    for i = 1 : numBetas
        alpha = alphaMat(i,j);

        beta = betaVec(i);

        polInd = 1;

        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static level pol';

        if i == 1
            h = scatter(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ1(:, colum_no),  resP1(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color);
        end

        plot(resEQ1(:, colum_no), resP1(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1.2);


        polInd = 2;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Static time pol.';
        if i == 1
            h = scatter(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ2(:, colum_no),  resP2(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end

        h2= plot(resEQ2(:, colum_no), resP2(:, colum_no)./ (1 - exp(-alpha /beta)), 'DisplayName', ['\beta= ' ,  num2str(beta)], 'Color', color,'LineWidth',1.2);
        ratioLegendEntries = [ratioLegendEntries, h2];
        ratioLegendLabels{end+1} =  ['\beta= ' ,  num2str(beta)];

        polInd = 3;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta-pol.';
        if i == 1
            h = scatter(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ3(:, colum_no),  resP3(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ3(:, colum_no), resP3(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);


        polInd = 4;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online c-pol.';
        if i == 1
            h = scatter(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ4(:, colum_no),  resP4(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ4(:, colum_no), resP4(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);

        polInd =5;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^2-pol.';
        if i == 1
            h = scatter(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ5(:, colum_no),  resP5(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ5(:, colum_no), resP5(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);
        


        polInd =6;
        baseColor = baseRatioColors{polInd};

        color = baseColor * shadeFactor(i);
        color = min(max(color, 0), 1); % Ensure color values are within valid range
        name= 'Online \theta\cdotq^{0.5}-pol.';
        if i == 1
            h = scatter(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'DisplayName', name, 'MarkerEdgeColor', color,'LineWidth',0.7);
            legendEntries = [legendEntries, h];
            legendLabels{end+1} = name;
        else
            scatter(resEQ6(:, colum_no),  resP6(:, colum_no)./ (1 - exp(-alpha /beta)), policyMarkers{polInd}, 'MarkerEdgeColor', color,'LineWidth',0.7);
        end
        plot(resEQ6(:, colum_no), resP6(:, colum_no)./ (1 - exp(-alpha /beta)), 'Color', color','LineWidth',1);


        colum_no=        colum_no+3;

    end
    hold off;
    legend([legendEntries, ratioLegendEntries], [legendLabels, ratioLegendLabels],'Location', 'best'); % Combine policy and ratio legends
    title(['Pareto Frontiers for \alpha/\beta = ', num2str(ratio)]);
    xlabel('Mean Number Waiting')
    ylabel('Premature Closure Rate for Active Customers')
    set(gca, 'yscale', 'log');
    set(gca, 'xscale', 'log');
    set(gca, 'FontSize', 14);
    name =  sprintf('V2_Stc_Onl6_ParetoF_rat%.1f_lmb%.1f.fig', ratio, lambda);
    name2 = sprintf('V2_Stc_Onl6_ParetoF_rat%.1f_lmb%.1f.jpg',  ratio, lambda);
    savefig(name);
    saveas(gcf, name2);
    close;
end








