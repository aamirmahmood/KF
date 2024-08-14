
clc; close all;clear all;


delta_T = 1; %measurements period
sigma_a = 0.2; %The random acceleration standard deviation (m/s)
sigma_r_m = 5; %The range measurement error standard deviation (m)
sigma_phi_m = 0.0087; %The bearing angle measurement error standard deviation (rad)


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


%% EKF Filtering


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



% Creating vehicle position plot

% Extract x and y estimates
x_estimate_values = x_estimates(1, :);
y_estimate_values = x_estimates(4, :);



% Absolute error - X coordinate
figure(1);
absError_xCoord_EKF =  abs(x_total - x_estimate_values);
plot([1:35], absError_xCoord_EKF, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
hold on;


% Absolute error - Y coordinate
figure(2);
absError_yCoord_EKF =  abs(y_total - y_estimate_values);
plot([1:35], absError_yCoord_EKF, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
hold on;


% Uncertainty - X coordinate
figure(3);
plot([1:35], x_uncertainty, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
hold on;

% Uncertainty - Y coordinate
figure(4);
plot([1:35], y_uncertainty, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
hold on;



%% UKF

% UKF parameters computation

N = 6; % The number of dimensions
SP = 2*N + 1; %The number of sigma points: 2N + 1 = 13
kappa = 3 - N; % For Gaussian distribution, the rule of thumb is to set: N + κ = 3
N_plus_kappa = 3;

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

    % 1. UKF-step 1: Calculate Sigma Points
        x_SP_zero = x_estimate_init; % x_SP^0: the first Sigma Point
        Cholesky_matrix = chol(N_plus_kappa*uncertainty_estimate_init); % or 3* diag ( repmat (500 ,1 ,6) 
        x_SP_init = [x_SP_zero, x_SP_zero + Cholesky_matrix, x_SP_zero - Cholesky_matrix];

    % 2. Points propagation
        % X(1,0) = f(X(0,0)). Since the dynamic model is linear, we can
        % write: X(1,0) = F * X(0,0)
       x_SP_pred = F*x_SP_init;

    % 3. Weights Calculation   
      omega_0 = kappa/N_plus_kappa;
      omega_1 = 1/(2*N_plus_kappa);
      omega = [omega_0 ones(1,2*N)*omega_1]';

    % 4. Mean and covariance computation  
      x_estimate_pred = x_SP_pred*omega;
      uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(omega)*(x_SP_pred - x_estimate_pred)' + Q; 

% Storage for estimates and uncertainties
x_estimates = zeros(6, length(r));
uncertainty_estimates = zeros(6, 6, length(r));
%x_y_uncertainty = zeros(2,2, length(r));

% KF Loop 
for ii = 1:length(r)
    
    % 1. Measure
    z = [r(ii); phi(ii)];

    % 2. Update
    % Using the Unscented Transform, transfer the state sigma points X(1,0) to the 
    % measurement space Z using the measurement function h(x):
    
    % In order to compute Z1, we perform elementwise operations between the first row
    % of X(1,0) that represents the x position (X(1,0) (1, :)) and the fourth row of X(1,0) that
    % represents the y position (X(1,0) (4, :)):

    Z = [sqrt(x_SP_pred(1,:).^2 + x_SP_pred(4,:).^2); atan(x_SP_pred(4,:)./x_SP_pred(1,:))];

    z_mean = Z * omega;

    % Compute covariance at the measurement space    
    
    Pz  = (Z-z_mean)*diag(omega)*(Z-z_mean)' + R;
    Pxz = (x_SP_pred - x_estimate_pred)*diag(omega)*(Z-z_mean)';

    % Kalman gain 
    K = Pxz/Pz;

    % Update estimate with measurement:
    x_estimate = x_estimate_pred + K*(z - z_mean);

    % Update the current uncertainty estimate (covariance of the estimate)
    uncertainty_estimate = uncertainty_estimate_pred - K*Pz*K';
    
    % Store the estimates and uncertainties
    x_estimates(:, ii) = x_estimate;
    uncertainty_estimates(:, :, ii) = uncertainty_estimate;
    %x_y_uncertainty(:, :, ii) = [uncertainty_estimate(1,1) 0; 0 uncertainty_estimate(4,4)];

    % 3. Predict 
        % 1. UKF-step 1: Calculate Sigma Points
        x_SP_zero = x_estimate; % x_SP^0: the first Sigma Point
        Cholesky_matrix = chol(N_plus_kappa*uncertainty_estimate)';
        x_SP_init = [x_SP_zero, x_SP_zero + Cholesky_matrix, x_SP_zero - Cholesky_matrix];

    % 2. Points propagation
        % X(1,0) = f(X(0,0)). Since the dynamic model is linear, we can
        % write: X(1,0) = F * X(0,0)
       x_SP_pred = F*x_SP_init;

    % 4. Mean and covariance computation  
      x_estimate_pred = x_SP_pred*omega;
      uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(omega)*(x_SP_pred - x_estimate_pred)' + Q; 


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


% Creating vehicle position plot

% Extract x and y estimates
x_estimate_values = x_estimates(1, :);
y_estimate_values = x_estimates(4, :);


% Absolute error - X coordinate
figure(1);
absError_xCoord_UKF =  abs(x_total - x_estimate_values);
plot([1:35], absError_xCoord_UKF, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel(' |X error|(m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Absolute error - X coordinate', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
legend('EKF', 'UKF')
% Customizing the font size for axes
set(gca, 'FontSize', 14);


% Absolute error - Y coordinate
figure(2);
absError_yCoord_UKF =  abs(y_total - y_estimate_values);
plot([1:35], absError_yCoord_UKF, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel(' |Y error|(m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Absolute error - Y coordinate', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
legend('EKF', 'UKF')
% Customizing the font size for axes
set(gca, 'FontSize', 14);


% Uncertainty - X coordinate
figure(3);
plot([1:35], x_uncertainty, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel(' \sigma_x (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Uncertainty - X coordinate', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
legend('EKF', 'UKF')
% Customizing the font size for axes
set(gca, 'FontSize', 14);

% Uncertainty - Y coordinate
figure(4);
plot([1:35], y_uncertainty, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel(' \sigma_y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Uncertainty - Y coordinate', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
legend('EKF', 'UKF')
% Customizing the font size for axes
set(gca, 'FontSize', 14);