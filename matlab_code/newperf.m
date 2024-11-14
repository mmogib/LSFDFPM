function succeed = newperf(A, mt, varargin)
    % Evaluate performance and plot curves with enhanced visuals
    % ibrahimkarym@gmail.com

    [n_p, n_s] = size(A);  % Get matrix dimensions

    % Use a color map for better visual differentiation
    cmap = lines(n_s);  % 'lines' provides distinct colors for multiple lines
    Ls = ["-", "--", "-.", ":", "-", "--", "-.", ":"];  % Line styles

    % Normalize performance matrix
    R = zeros(n_p, n_s);
    for i = 1:n_p
        R(i, :) = A(i, :) / min(A(i, :));
    end

    % Calculate success rate
    succeed = sum(R == 1) / n_p;

    % Plot setup
    figure;
    hold on;
    tou = 1:0.1:mt;  % Range of tolerance values
    rho = zeros(n_s, numel(tou));  % Preallocate rho for performance ratios

    % Generate performance curves for each solver
    for s = 1:n_s
        for j = 1:numel(tou)
            T = tou(j) * ones(n_p, 1);  % Create threshold vector
            rho(s, j) = sum(R(:, s) <= T) / n_p;  % Calculate ratio
        end
        % Plot each curve with distinct colors and styles
        plot(tou, rho(s, :), 'Color', cmap(s, :), ...
            'LineStyle', Ls(mod(s - 1, numel(Ls)) + 1), ...
            'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
            'DisplayName', varargin{s});
    end

    % Configure plot appearance
    plt_title = 'Performance Profiles';
    if length(varargin)> n_s
        plt_title = varargin{end};
    end
    legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'none');
    axis([1 mt 0 1]);  % Set axis limits
    grid on;  % Add grid
    grid minor;  % Add minor grid lines
    xlabel('\rho', 'FontSize', 14);  % X-axis label
    ylabel('Performance Ratio', 'FontSize', 14);  % Y-axis label
    title(plt_title, 'FontSize', 16);  % Title
    hold off;
end
