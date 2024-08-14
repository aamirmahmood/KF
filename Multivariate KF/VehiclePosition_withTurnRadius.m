% Define the initial position
x_initial = 300;
y_initial = -400;

% Define the final position after the turn
x_final = 20;
y_final = 300;

% Define the number of points
num_points_straight = 16; % Number of points for the straight line
num_points_turn = 19; % Number of points for the turn

% Define the straight line segment
x_straight = x_initial * ones(1, num_points_straight);
y_straight = linspace(y_initial, 0, num_points_straight);

% Calculate the angle required for the turn
turn_radius = 300;
arc_length = sqrt((x_initial - x_final)^2 + (y_final - 0)^2);
theta_final = 2 * asin(arc_length / (2 * turn_radius)); % Calculate theta required to reach final position

% Define the turning maneuver (left turn with 300 meters radius)
theta = linspace(0, theta_final, num_points_turn); % Angle from 0 to theta_final radians
x_turn = x_initial - turn_radius * (1 - cos(theta)); % X coordinates during the turn
y_turn = turn_radius * sin(theta); % Y coordinates during the turn

% Combine straight and turning parts
x_total = [x_straight, x_turn];
y_total = [y_straight, y_turn];

% Plot the trajectory
figure;
plot(x_total, y_total, 'g-', 'LineWidth', 2);
hold on;
plot(x_total, y_total, 'bo-', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle Position', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)

% Set figure size to be similar to the uploaded image
set(gcf, 'Position', [100, 100, 400, 800]);

% Ensure equal scaling for both axes
axis equal;

% Customizing the font size for axes
set(gca, 'FontSize', 14);

% Save the figure as a PNG file
saveas(gcf, 'vehicle_position.png');
