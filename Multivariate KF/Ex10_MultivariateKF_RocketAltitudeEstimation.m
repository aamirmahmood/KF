
clc; close all;clear all;

%% Model

% Parameters
dt = 0.25; % measurement period in seconds
acc = 30; % acceleration in m/s^2
total_time = 7.25; % total time in seconds (adjust as needed)
num_points = total_time / dt; % number of data points

% Initialize vectors
t = linspace(0, total_time, num_points); % time vector
h_true = zeros(1, num_points); % height vector
v_true = zeros(1, num_points); % height vector

% Initial conditions
h0 = 15; % initial height in meters
v0 = 0; % initial velocity in m/s

% Calculate height at each time step
for i = 1:num_points
    h_true(i) = h0 + v0 * t(i) + 0.5 * acc * t(i)^2;
end

% Calculate true velocity at each time step
for i = 1:num_points
    v_true(i) = v0 + acc * t(i);
end

%% Data

% Define the data
% Define the vectors h and a
h = [6.43, 1.3, 39.43, 45.89, 41.44, 48.7, 78.06, 80.08, ...
     61.77, 75.15, 110.39, 127.83, 158.75, 156.55, 213.32, 229.82, 262.8, ...
     297.57, 335.69, 367.92, 377.19, 411.18, 460.7, 468.39, 553.9, 583.97, ...
     655.15, 723.09, 736.85, 787.22];

a = [39.81, 39.67, 39.81, 39.84, 40.05, 39.85, 39.78, 39.65, ...
     39.67, 39.78, 39.59, 39.87, 39.85, 39.59, 39.84, 39.9, 39.63, ...
     39.59, 39.76, 39.79, 39.73, 39.93, 39.83, 39.85, 39.94, 39.86, ...
     39.76, 39.86, 39.74, 39.94];


%% Filtering

delta_T = 0.25; %measurements period (s)
rocket_a = 30; %The random acceleration (m/s^2)
sigma_x_m = 20; %The altimeter measurement error standard deviation (m)
epsilon = 0.1; %The accelerometer measurement error standard deviation (m/s)

% Transition matrix (F)
F = [1 delta_T;...
     0 1      ];

G = [0.5*delta_T^2;...
     delta_T];

u = 0;


% Process noise matrix (Q)
Q = epsilon^2*[delta_T^4/4 delta_T^3/2;...
               delta_T^3/2 delta_T^2  ];

% Measurment uncertainty 
R = [sigma_x_m^2];


% Iteration Zero
% 1. Initialization 

x_estimate_init = [0 0]'; %x^hat(0,0)
uncertainty_estimate_init = [500 0  ;...
                             0   500]; % uncertainty estimate P(0,0)
% 2. Predication
x_estimate_pred           = F*x_estimate_init + G*u;
uncertainty_estimate_pred = F*uncertainty_estimate_init*F' + Q; 

% Storage for estimates and uncertainties
x_estimates = zeros(length(x_estimate_init), length(h));
uncertainty_estimates = zeros(length(x_estimate_init), length(x_estimate_init), length(h));
%x_y_uncertainty = zeros(2,2, length(x));

% KF Loop 
for ii = 1:length(h)
    % 1. Measure
    z = [h(ii)];
    u = [a(ii)];
    %The dimension of zn is 2×1 and the dimension of xn is 6×1. Therefore, the dimension 
    % of the observation matrix H shall be 2 × 6.
    H = [1 0];

    % 2. Update
    % Kalman gain 
    K = uncertainty_estimate_pred * H' * (H * uncertainty_estimate_pred * H' + R)^(-1); 

    % Estimate the current state
    x_estimate = x_estimate_pred + K*(z - H*x_estimate_pred);

    % Update the current uncertainty estimate
    uncertainty_estimate = (eye(2) - K*H)*uncertainty_estimate_pred * (eye(2) - K*H)' + K * R * K';
    
    % Store the estimates and uncertainties
    x_estimates(:,ii) = x_estimate;
    uncertainty_estimates(:, :, ii) = uncertainty_estimate;
    %x_y_uncertainty(:, :, ii) = [uncertainty_estimate(1,1) 0; 0 uncertainty_estimate(4,4)];

    % 3. Predict 
    x_estimate_pred = F*x_estimate + G* (u + (-9.8));
    uncertainty_estimate_pred = F * uncertainty_estimate * F' + Q;

end

%% Confidence ellipse

CI = 95;
temp = (1 - CI/100)/2;

uncertainty_x = squeeze(uncertainty_estimates(1,1,:))';
uncertainty_y = squeeze(uncertainty_estimates(2,2,:))';


x_low = x_estimates(1,:) - norminv(temp)*sqrt(uncertainty_x);
x_up  = x_estimates(1,:) + norminv(temp)*sqrt(uncertainty_x);

y_low = x_estimates(2,:) - norminv(temp)*sqrt(uncertainty_y);
y_up  = x_estimates(2,:) + norminv(temp)*sqrt(uncertainty_y);


%% Plotting


% Creating rocket altitude plot


figure;
plot(t, h_true, 'gs-', 'LineWidth', 2); % Measurements 
x_axis = [0:delta_T:(length(h)-1)*delta_T];
hold on;
plot(x_axis, h, 'bs-', 'LineWidth', 2); % Measurements 
plot(x_axis, x_estimates(1,:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Altitude (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Rocket altitude', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)


plot(x_axis, x_low, 'r', 'LineWidth',1);
plot(x_axis, x_up, 'r', 'LineWidth',1);
x2 = [x_axis, fliplr(x_axis)];
inBetween = [x_low, fliplr(x_up)];
p4 = patch(x2, inBetween, 'y');
%legend(p1,'Area 1')
alpha(0.25)

legend('True values', 'Measurements', 'Estimates', '95% CI')


% Customizing the font size for axes
set(gca, 'FontSize', 14);

%2. Creating rocket x-axis/y-axis velocity plot

figure;
plot(t, v_true, 'gs-', 'LineWidth', 2); % Measurements 
x_axis = [0:delta_T:(length(h)-1)*delta_T];
hold on;
plot(x_axis, x_estimates(2,:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Velocity (m/s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Rocket Velocity', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)


plot(x_axis, y_low, 'r', 'LineWidth',1);
plot(x_axis, y_up, 'r', 'LineWidth',1);
x2 = [x_axis, fliplr(x_axis)];
inBetween = [y_low, fliplr(y_up)];
p4 = patch(x2, inBetween, 'y');
%legend(p1,'Area 1')
alpha(0.25)

legend('True values', 'Estimates', '95% CI')


% Customizing the font size for axes
set(gca, 'FontSize', 14);

