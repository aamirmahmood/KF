
clc; close all;clear all;


%% Model
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
y_straight = linspace(y_initial, 0.01, num_points_straight);

% Calculate the angle required for the turn
turn_radius = 300;
arc_length = sqrt((x_initial - x_final)^2 + (y_final - 0.01)^2);
theta_final = 2 * asin(arc_length / (2 * turn_radius)); % Calculate theta required to reach final position

%TODO based on new y_init = 0.01, calculate the right starting theta

% Define the turning maneuver (left turn with 300 meters radius)
theta = linspace(0.085, theta_final, num_points_turn); % Angle from 0 to theta_final radians
x_turn = x_initial - turn_radius * (1 - cos(theta)); % X coordinates during the turn
y_turn = turn_radius * sin(theta); % Y coordinates during the turn

% Combine straight and turning parts
x_total = [x_straight, x_turn];
y_total = [y_straight, y_turn];

%% Data

% Define the data
r = [502.55, 477.34, 457.21, 442.94, 427.27, 406.05, 400.73, 377.32, ...
     360.27, 345.93, 333.34, 328.07, 315.48, 301.41, 302.87, 304.25, ...
     294.46, 294.29, 299.38, 299.37, 300.68, 304.10, 301.96, 300.30, ...
     301.90, 296.70, 297.07, 295.29, 296.31, 300.62, 292.30, 298.11, ...
     298.07, 298.92, 298.04];
phi = [-0.9316, -0.8977, -0.8512, -0.8114, -0.7853, -0.7392, -0.7052, -0.6478, ...
       -0.5900, -0.5183, -0.4698, -0.3952, -0.3026, -0.2445, -0.1626, -0.0937, ...
        0.0085,  0.0856,  0.1675,  0.2467,  0.3290,  0.4149,  0.5040,  0.5934, ...
        0.6670,  0.7537,  0.8354,  0.9195,  1.0039,  1.0923,  1.1546,  1.2564, ...
        1.3274, 1.409, 1.5011];


%% Filtering

delta_T = 1; %measurements period
sigma_a = 0.2; %The random acceleration standard deviation (m/s)
sigma_r_m = 5; %The range measurement error standard deviation (m)
sigma_phi_m = 0.0087; %The bearing angle measurement error standard deviation (rad)

% Transition matrix (F)
F = [1 delta_T  0.5*delta_T   0         0       0;...
     0 1        delta_T       0         0       0;...
     0 0        1             0         0       0;...
     0 0        0             1         delta_T 0.5*delta_T;...
     0 0        0             0         1       delta_T;...
     0 0        0             0         0       1]; % uncertainty estimate P(0,0)

% Process noise matrix (Q)
Q = sigma_a^2*[delta_T^4/4 delta_T^3/2 delta_T^2/2   0           0            0;...
     delta_T^3/2 delta_T^2   delta_T       0           0            0;...
     delta_T^2/2 delta_T     1             0           0            0;...
     0           0           0             delta_T^4/4 delta_T^3/2  delta_T^2/2;...
     0           0           0             delta_T^3/2 delta_T^2    delta_T;...
     0           0           0             delta_T^2/2 delta_T      1];

% Measurment uncertainty 
R = [sigma_r_m^2 0; 
     0           sigma_phi_m^2];

% Define the data
r = [502.55, 477.34, 457.21, 442.94, 427.27, 406.05, 400.73, 377.32, ...
     360.27, 345.93, 333.34, 328.07, 315.48, 301.41, 302.87, 304.25, ...
     294.46, 294.29, 299.38, 299.37, 300.68, 304.10, 301.96, 300.30, ...
     301.90, 296.70, 297.07, 295.29, 296.31, 300.62, 292.30, 298.11, ...
     298.07, 298.92, 298.04];
phi = [-0.9316, -0.8977, -0.8512, -0.8114, -0.7853, -0.7392, -0.7052, -0.6478, ...
       -0.5900, -0.5183, -0.4698, -0.3952, -0.3026, -0.2445, -0.1626, -0.0937, ...
        0.0085,  0.0856,  0.1675,  0.2467,  0.3290,  0.4149,  0.5040,  0.5934, ...
        0.6670,  0.7537,  0.8354,  0.9195,  1.0039,  1.0923,  1.1546,  1.2564, ...
        1.3274, 1.409, 1.5011];

% Iteration Zero
% 1. Initialization 

%We don%t know the vehicle location, so we approximate the initial position at about 
% 100m from the true vehicle position (ˆx0,0 = 400m, ˆy0,0 = −300m)

x_estimate_init = [400 0 0 -300 0 0]'; %x^hat(0,0)
uncertainty_estimate_init = [500 0   0   0    0   0;...
                             0   500 0   0    0   0;...
                             0   0   500 0    0   0;...
                             0   0   0   500  0   0;...
                             0   0   0   0    500 0;...
                             0   0   0   0    0   500]; % uncertainty estimate P(0,0)
% 2. Predication
x_estimate_pred           = F*x_estimate_init;
uncertainty_estimate_pred = F*uncertainty_estimate_init*F' + Q; 

% Storage for estimates and uncertainties
x_estimates = zeros(6, length(r));
uncertainty_estimates = zeros(6, 6, length(r));
x_y_uncertainty = zeros(2,2, length(r));

% KF Loop 
for ii = 1:length(r)
    
    % 1. Measure
    z = [r(ii); phi(ii)];

    % 2. Update

    %Observation matrix (h(ˆx1,0)) calculation.
    x = x_estimate_pred(1,1);
    y = x_estimate_pred(4,1);

    h_x          = [sqrt(x^2 + y^2); atan(y/x)];

    h_x_jacobian =  [x/sqrt(x^2 + y^2) 0 0 y/sqrt(x^2 + y^2) 0 0; ...
                     -y/(x^2 + y^2) 0 0 x/(x^2 + y^2) 0 0];

    % Kalman gain 
    K = uncertainty_estimate_pred * h_x_jacobian' * (h_x_jacobian * uncertainty_estimate_pred * h_x_jacobian' + R)^(-1); 

    % Estimate the current state
    x_estimate = x_estimate_pred + K*(z - h_x);

    % Update the current uncertainty estimate
    uncertainty_estimate = (eye(6) - K*h_x_jacobian)*uncertainty_estimate_pred * (eye(6) - K*h_x_jacobian)' + K * R * K';
    
    % Store the estimates and uncertainties
    x_estimates(:, ii) = x_estimate;
    uncertainty_estimates(:, :, ii) = uncertainty_estimate;
    %x_y_uncertainty(:, :, ii) = [uncertainty_estimate(1,1) 0; 0 uncertainty_estimate(4,4)];

    % 3. Predict 
    x_estimate_pred = F*x_estimate;
    uncertainty_estimate_pred = F * uncertainty_estimate * F' + Q;

end

%% Confidence ellipse
CI = 95; % in percentage 

x_y_estimate(:,:) = [x_estimates(1,:);x_estimates(4,:)];

% Extract the uncertainty for x (1,1 element) and y (4,4 element)
x_uncertainty = squeeze(uncertainty_estimates(1, 1, :)); % Uncertainty in x
y_uncertainty = squeeze(uncertainty_estimates(4, 4, :)); % Uncertainty in y
% Combine x and y uncertainties into a 2x2xN array
x_y_uncertainty = zeros(2, 2, length(r));
x_y_uncertainty(1, 1, :) = x_uncertainty;
x_y_uncertainty(2, 2, :) = y_uncertainty;

r_ellipses = calculate_error_ellipses (x_y_estimate, x_y_uncertainty, CI);


%% Plotting

%% estimation uncertainty 
figure(20);

plot([1:35], x_uncertainty, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
hold on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Varinace (m^2)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%title('Vehicle x-axis velocity', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
%legend('True values', 'Estimates')

% Creating vehicle position plot

% Extract x and y estimates
x_estimate_values = x_estimates(1, :);
y_estimate_values = x_estimates(4, :);

figure;
plot(x_total, y_total, 'gd-', 'LineWidth', 1.5,'Color', [0.1, 0.63, 0.1]); % Model based--true values
hold on;
plot(r.*cos(phi), r.*sin(phi), 'bs-', 'LineWidth', 2); % Measurements 
plot(x_estimate_values, y_estimate_values, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle Position', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
legend('True values', 'Measurements', 'Estimates')

% Set figure size
set(gcf, 'Position', [100, 100, 400, 800]);

% Ensure equal scaling for both axes
axis equal;

% Customizing the font size for axes
set(gca, 'FontSize', 14);

%2. Creating vehicle x-axis/y-axis velocity plot

% Extract x- and y- velocity estimates
v_x_estimate = x_estimates(2, :);
v_y_estimate = x_estimates(5, :);



%True Velocity
v_x_true = diff(x_total)/delta_T;
v_y_true = diff(y_total)/delta_T;


figure;
plot([1:34], v_x_true, 'gd-', 'LineWidth', 1.5,'Color', [0.1, 0.63, 0.1]); % Model based--true values
hold on;
plot([1:34], v_x_estimate(1:34), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle x-axis velocity', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
legend('True values', 'Estimates')


figure;
plot([1:34], v_y_true, 'gd-', 'LineWidth', 1.5,'Color', [0.1, 0.63, 0.1]); % Model based--true values
hold on;
plot([1:34], v_y_estimate(1:34), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle y-axis velocity', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
legend('True values', 'Estimates')



% plot with error ellipses


% Creating vehicle position plot

% Extract x and y estimates
x_estimate_values = x_estimates(1, :);
y_estimate_values = x_estimates(4, :);

figure;
plot(x_total, y_total, 'gd-', 'LineWidth', 1.5,'Color', [0.1, 0.63, 0.1]); % Model based--true values
hold on;
plot(r.*cos(phi), r.*sin(phi), 'bs-', 'LineWidth', 2); % Measurements 
plot(x_estimate_values, y_estimate_values, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate

for ii = 1:length(r_ellipses)
    r_ellipse = r_ellipses{ii};
    fill(r_ellipse(:, 1), r_ellipse(:, 2), 'yellow', 'EdgeColor', [0.8, 0.8, 0], 'FaceAlpha', 0.5);

    %plot(r_ellipse(:, 1), r_ellipse(:, 2), '-');
end

grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle Position', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
legend('True values', 'Measurements', 'Estimates', '95% Confidence Interval')

% Set figure size
set(gcf, 'Position', [100, 100, 400, 800]);

% Ensure equal scaling for both axes
axis equal;

% Customizing the font size for axes
set(gca, 'FontSize', 14);



%% 
% plot with error ellipses


% Creating vehicle position plot

% Extract x and y estimates
x_estimate_values = x_estimates(1, :);
y_estimate_values = x_estimates(4, :);

figure;
plot(x_total, y_total, 'gd-', 'LineWidth', 1.5,'Color', [0.1, 0.63, 0.1]); % Model based--true values
hold on;
plot(r.*cos(phi), r.*sin(phi), 'bs-', 'LineWidth', 2); % Measurements 
plot(x_estimate_values, y_estimate_values, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate

for ii = 1:length(r_ellipses)
    r_ellipse = r_ellipses{ii};
    fill(r_ellipse(:, 1), r_ellipse(:, 2), 'yellow', 'EdgeColor', [0.8, 0.8, 0], 'FaceAlpha', 0.5);

    %plot(r_ellipse(:, 1), r_ellipse(:, 2), '-');
end

grid on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Vehicle Position', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)
legend('True values', 'Measurements', 'Estimates', '95% Confidence Interval')

% Set figure size
set(gcf, 'Position', [100, 100, 400, 800]);

% Ensure equal scaling for both axes
axis equal;

% Customizing the font size for axes
set(gca, 'FontSize', 14);
axis([230 320 -400 100])


