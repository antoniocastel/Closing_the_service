
%% Code to generate Figure 1 of CLOSING THE SERVICE: Contrasting Activity-Based and Time-Based Systematic Closure Policies
% By Antonio Castellanos, Andrew Daw, Amy Ward and Galit B. Yom-Tov

    % Parameters for the Hawkes process
    alpha = 7.9/60;   % influence of past events
    beta = 8.15/60;    % decay of intensity
    

    T = 500;       % total time to simulate

    % Initialize with the first event at time 0
    events = [0];  
    t = 0;         % Start from time zero

    while t < T
        % Calculate the conditional intensity function \mu_t
        mu_t = sum(alpha * exp(-beta * (t - events)));

        % Generate the next event time with the thinning algorithm
        u = rand();  % Uniformly distributed random number 
        tau = -log(u) / mu_t; % Time until next event
        t = t + tau; % Proposed next event time
        
        % Ensure we do not go past the end time T
        if t >= T
            break; % If the proposed time is outside the simulation window, break
        end

        % Determine if the event is accepted
        mu_t_prime = sum(alpha * exp(-beta * (t - events)));
        d = rand();  % uniformly distributed random number
        if d <= mu_t_prime / mu_t
            events = [events t]; % Accept  event
        end
    end

    % Calculate intensity  for plotting
    time_points = linspace(0, T, 1000000);
    intensity_values = arrayfun(@(t) sum(alpha * exp(-beta * (t - events(events < t)))), time_points);

    % Plot the intensity and the event arrivals
    figure;
    plot(time_points, intensity_values, 'b-', 'LineWidth', 1.5); hold on;
    scatter(events, zeros(size(events)), 50, 'or', 'filled'); % Show arrivals as red dots
    xlim([0 120])
    xlabel('Time');
    ylabel('Intensity');
    xticks(0:10:T);
    set(gca,'FontSize',14)
    yline(0.04, 'k--', 'LineWidth', 1.5);
        xline(25, '--', 'LineWidth', 1.5, 'Color',[0.0000, 0.4275, 0.1725]);
    yLimits = get(gca, 'YLim'); 
    fill([15.3 15.3 25 25], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.0000, 0.4275, 0.1725], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
   
    legend('\mu_t', 'Events', 'Activity level policy', 'Inactive time policy','Location', 'northeast');
    text(52, 0.052, '$l$', 'FontSize', 17, 'HorizontalAlignment', 'right', Interpreter='latex', FontWeight='bold'); 
    text(52, 0.052, '$l$', 'FontSize', 17, 'HorizontalAlignment', 'right', Interpreter='latex', FontWeight='bold');
    text(21, 0.46, '\epsilon', 'FontSize', 17, 'HorizontalAlignment', 'right'); % Adjust the -5 according to your figure's actual margins

    hold off;
