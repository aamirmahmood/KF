clc; close all; clear all;

UKF_modified = 1; % Set it to 0 when using original UKF, 1 if modified sigma point in UKF


%% Model
% Define the parameters
L = 0.5; % Length of the pendulum (m)
g = 9.8; % Gravitational acceleration (m/s^2)
theta0 = 0.2; % Initial angle (radians)
omega0 = 0; % Initial angular velocity (rad/s)
sigma_a = 1; % Process Noise Uncertainty (rad/s^2)
T = 2.5; % Total time (s)
dt = 0.05; % Time step (s)
t = 0:dt:T; % Time vector

% Initialize arrays
theta = zeros(size(t));
omega = zeros(size(t));

% Set initial conditions
theta(1) = theta0;
omega(1) = omega0;

% Generate process noise for angular acceleration
a_noise = sigma_a * randn(size(t));

% Simulate the pendulum with process noise
for i = 2:length(t)
    % Update angular velocity with process noise
    omega(i) = omega(i-1) - (g/L) * sin(theta(i-1)) * dt + a_noise(i) * dt;
    % Update angle
    theta(i) = theta(i-1) + omega(i) * dt;
end


% Plot theta over time
figure;
subplot(3,1,1);
plot(t, theta, 'r-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('\theta (rad)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Pendulum Angle', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;

% Plot angular velocity over time
subplot(3,1,2);
plot(t, omega, 'g-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Angular Velocity (rad/s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Pendulum Angular Velocity', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;

% Plot x over time
% Compute x = L * sin(theta)
x_pos = L * sin(theta);
subplot(3,1,3);
plot(t, x_pos, 'b-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('x (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Pendulum Displacement', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;

%% Data

% Define the data, adding some random noise to the measurement as full
% measurements are not given in the book and we need to create x from the
% true theta generated above
%x = [0.119 0.113 0.12 0.101 0.099 0.063 0.008 -0.017 -0.037 -0.05]; % data in the book
x = x_pos + 0.00001*randn(1,length(x_pos));



%% Filtering

delta_T = dt; % measurements period
L = 0.5; %The Pendulum string length (m)
g = 9.8; % Gravitational acceleration constant (m/s2)
sigma_x_m = 0.01; % Measurement uncertainty (standard deviation) (m)
sigma_a = 1; % Process Noise Uncertainty (angular acceleration standard deviation) (rad/s)

% Process noise matrix (Q)
Q = sigma_a^2*[delta_T^4/4 delta_T^3/2;...
               delta_T^3/2 delta_T^2];

% Measurment uncertainty 
R = [sigma_x_m^2];

% UKF parameters computation

N = 2; % The number of dimensions
SP = 2*N + 1; %The number of sigma points: 2N + 1 = 5

% Based on modified sigma point computation algorithm, set:
kappa = 0;
alpha = 0.1;
beta  = 2;

lambda = alpha^2 * (N + kappa) - N;
N_plus_lambda = N + lambda;

if (UKF_modified==0)
    kappa = 3 - N; % For Gaussian distribution, the rule of thumb is to set: N + κ = 3
    N_plus_kappa = 3;
end



% Iteration Zero
% 1. Initialization 

x_estimate_init = [0.0873 0]'; %x^hat(0,0) % We don t know the pendulum angle, so our initial angle approximation includes an error. We set the initial velocity to 0.
uncertainty_estimate_init = [5 0; ...
                             0 5]; % uncertainty estimate P(0,0) % Since our initial state vector is a guess, we set a high estimate uncertainty. The high estimate uncertainty results in a high Kalman Gain by giving high weight to the measurement.

% 2. Predication

    % 1. UKF-step 1: Calculate Sigma Points
        x_SP_zero = x_estimate_init; % x_SP^0: the first Sigma Point
        Cholesky_matrix = chol(N_plus_lambda*uncertainty_estimate_init); % or 3* diag ( repmat (500 ,1 ,6) 
        if (UKF_modified==0)
           Cholesky_matrix = chol(N_plus_kappa*uncertainty_estimate_init)'; % or 3* diag ( repmat (500 ,1 ,6) 
        end
        x_SP_init = [x_SP_zero, x_SP_zero + Cholesky_matrix, x_SP_zero - Cholesky_matrix];

    % 2. Points propagation      
        x_SP_pred = [x_SP_init(1,:) + x_SP_init(2,:)*delta_T; ...
                     x_SP_init(2,:) - g/L * sin(x_SP_init(1,:))*delta_T]; % earlier in ex11, I called it x_estimate_pred   

    % 3. Weights Calculation   
        omega_0_m = lambda/(N+lambda);
        omega_0_c = lambda/(N+lambda) + (1 - alpha^2 + beta);
        omega_1   = 1/(2*(N+lambda));

        W_m = [omega_0_m ones(1,2*N)*omega_1];
        W_c = ([omega_0_c ones(1,2*N)*omega_1]);

        if (UKF_modified==0)
            omega_0 = kappa/N_plus_kappa;
            omega_1 = 1/(2*N_plus_kappa);
            omega_ukf = [omega_0 ones(1,2*N)*omega_1]';
        end
    % 4. Mean and covariance computation
      x_estimate_pred = x_SP_pred*W_m';
      uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(W_c)*(x_SP_pred - x_estimate_pred)' + Q; 
      if (UKF_modified==0)
          x_estimate_pred = x_SP_pred*omega_ukf;
          uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(omega_ukf)*(x_SP_pred - x_estimate_pred)' + Q;
      end


% Storage for estimates and uncertainties
x_estimates = zeros(2, length(x));
uncertainty_estimates = zeros(2, 2, length(x));
%x_y_uncertainty = zeros(2,2, length(r));



% KF Loop 
for ii = 1:length(x)
    
    % 1. Measure
    z = [x(ii)];

    % 2. Update
    % Using the Unscented Transform, transfer the state sigma points X(1,0) to the 
    % measurement space Z using the measurement function h(x):
    
    % In order to compute Z1, we perform elementwise operations between the first row
    % of X(1,0) that represents that represents the angle θ:

    Z = L*sin(x_SP_pred(1,:));

    z_mean = Z * W_m';

    if (UKF_modified==0)
        z_mean = Z * omega_ukf;
    end

    % Compute covariance at the measurement space    
    
    Pz  = (Z-z_mean)*diag(W_c)*(Z-z_mean)' + R;
    Pxz = (x_SP_pred - x_estimate_pred)*diag(W_c)*(Z-z_mean)';
    
    if (UKF_modified==0)
        Pz  = (Z-z_mean)*diag(omega_ukf)*(Z-z_mean)' + R;
        Pxz = (x_SP_pred - x_estimate_pred)*diag(omega_ukf)*(Z-z_mean)';
    end
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
        Cholesky_matrix = chol(N_plus_lambda*uncertainty_estimate)';
        if (UKF_modified==0)
            Cholesky_matrix = chol(N_plus_kappa*uncertainty_estimate)';
        end
        x_SP_init = [x_SP_zero, x_SP_zero + Cholesky_matrix, x_SP_zero - Cholesky_matrix];

    % 2. Points propagation
        x_SP_pred = [x_SP_init(1,:) + x_SP_init(2,:)*delta_T; ...
                     x_SP_init(2,:) - g/L * sin(x_SP_init(1,:))*delta_T]; % earlier in ex11, I called it x_estimate_pred

    % 4. Mean and covariance computation  
      x_estimate_pred = x_SP_pred*W_m';
      uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(W_c)*(x_SP_pred - x_estimate_pred)' + Q; 
      if(UKF_modified==0)
          x_estimate_pred = x_SP_pred*omega_ukf;
          uncertainty_estimate_pred = (x_SP_pred - x_estimate_pred)*diag(omega_ukf)*(x_SP_pred - x_estimate_pred)' + Q;
      end


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


% Creating pednulum angle plot
figure;
plot(t, theta, 'gs-', 'LineWidth', 2); % True angle values
hold on;
plot(t, asin(x/L), 'bs-', 'LineWidth', 2); % Measurements x = L*sin(theta), note that I add some randomness to x 
plot(t, x_estimates(1,:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Theta (rad)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Pendulum Angle', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)


plot(t, x_low, 'r', 'LineWidth',1);
plot(t, x_up, 'r', 'LineWidth',1);
x2 = [t, fliplr(t)];
inBetween = [x_low, fliplr(x_up)];
p4 = patch(x2, inBetween, 'y', 'FaceAlpha', 0.25);
%legend(p1,'Area 1')
%alpha(p4, 0.25);

legend('True values', 'Measurements', 'Estimates', '95% CI')


% Customizing the font size for axes
set(gca, 'FontSize', 14);

%2. Creating pendulum angular velocity plot

figure;
plot(t, omega, 'gs-', 'LineWidth', 2); % True angular values
hold on;
%plot(t, -sqrt(g/L) * theta0*sin(sqrt(g/L) * t), 'bs-', 'LineWidth', 2); % Measurements
plot(t, x_estimates(2,:), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6); % KF estimate
grid on;
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Velocity (m/s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
title('Pendulum Velocity', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
%set(gca, 'YDir', 'normal'); % Ensure the y-axis direction is normal (default)


plot(t, y_low, 'r', 'LineWidth',1);
plot(t, y_up, 'r', 'LineWidth',1);
x2 = [t, fliplr(t)];
inBetween = [y_low, fliplr(y_up)];
p4 = patch(x2, inBetween, 'y', 'FaceAlpha', 0.25);
%legend(p1,'Area 1')
%alpha(p4, 0.25)

%legend('True values', 'Measurements', 'Estimates', '95% CI'); % If Measurements plotting is enabled
legend('True values', 'Estimates', '95% CI')

% Customizing the font size for axes
set(gca, 'FontSize', 14);




