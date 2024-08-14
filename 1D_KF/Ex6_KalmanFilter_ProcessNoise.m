clc, close all, clear all;

T_true = 50;
z_true = [49.979 50.025 50 50.003 49.994 50.002 49.999 50.006 49.998 49.991];
z_true = [50.005 49.994 49.993 50.001 50.006 49.998 50.021 50.005 50 49.997]; % new values

% Measurements vector
z = [49.95 49.967 50.1 50.106 49.992 49.819 49.933 50.007 50.023 49.99];
z = [49.986 49.963 50.09 50.001 50.018 50.05 49.938 49.858 49.965 50.114];

% Measurement uncertainty 
sigma_std = 0.1; 
r = sigma_std^2;

% process noise variance
q = 0.0001;

%% ITERATION ZERO

% Initialization 

T_estimate = 10;                % Estimated Temperature
uncertainity_estimate = 10000;    % Estimate uncertainty with std. σ = 100

% Prediction at Iteration Zero
T_estimate_pred(1,1) = T_estimate;   % Since Dynamic Model is constant
uncertainity_estimate_pred(1,1) = uncertainity_estimate + q; % The extrapolated estimate uncertainty (variance) also doesn’t change

for ii = 1:length(z)
    
%   1) Measure
    z_n = z(ii); % New Measurement
    r_n = r;     % Measurement uncertainty
    
    
    % Kalman gain
    K(ii) = uncertainity_estimate_pred(1,ii) / (uncertainity_estimate_pred(1,ii) + r_n);
    
    
%   3) State Update
    T_estimate_current(ii)    = T_estimate_pred(1,ii) + K(ii)*(z(ii) - T_estimate_pred(1,ii)); %Estimating the current state
    uncertainity_estimate_current(ii) = (1-K(ii))* uncertainity_estimate_pred(1,ii); %Current estimate uncertainty
    
    % Prediction Updates
    T_estimate_pred(1, ii+1) = T_estimate_current (ii);
    uncertainity_estimate_pred(1,ii+1) = uncertainity_estimate_current(ii) + q;
end

x = 1:length(z);

figure();

subplot(3,1,1);
plot(x, z_true, '-d', 'LineWidth', 2);
hold on;
plot(x, z, '-*', 'LineWidth', 2);

plot(x, T_estimate_current, '-s', 'LineWidth', 2);
% plot([0], T_estimate, 'd', 'LineWidth', 4);

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Temperature (C)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('True Values', 'Measurements', 'Estimate', 'FontSize', 12, 'Location', 'NorthEast');
title('Temperature', 'FontSize', 12, 'FontWeight', 'Bold');


subplot(3,1,2);
% plot(x, ones(1,length(z))*r, '-d', 'LineWidth', 2);
% hold on;
plot(x, uncertainity_estimate_current, '-*', 'LineWidth', 2);

% plot(x, h_estimate_current, '-s', 'LineWidth', 2);
% plot([0], uncertainity_estimate, 'd', 'LineWidth', 4);

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Uncertainty', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('Estimate Uncertainty', 'FontSize', 12, 'Location', 'NorthEast');
title('Uncertainty', 'FontSize', 12, 'FontWeight', 'Bold');

subplot(3,1,3);

plot(x, K, '-d', 'LineWidth', 2);
xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
title('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');


%%
figure;

p1 = plot(x, z_true, '-d', 'LineWidth', 2);
hold on;
p2 = plot(x, z, '-*', 'LineWidth', 2);
p3 = plot(x, T_estimate_current, '-s', 'LineWidth', 2);
% 95% Confidence Interval

CI = 95;
temp = (1 - CI/100)/2;

low = T_estimate_current - norminv(temp)*sqrt(uncertainity_estimate_current);
up  = T_estimate_current + norminv(temp)*sqrt(uncertainity_estimate_current);

plot(x, low, 'r', 'LineWidth',1);
plot(x, up, 'r', 'LineWidth',1);
x2 = [x, fliplr(x)];
inBetween = [low, fliplr(up)];
p4 = patch(x2, inBetween, 'y');
%legend(p1,'Area 1')
alpha(0.25)
grid on;

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Height (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
%legend('True Values', 'Estimates', 'FontSize', 12, 'Location', 'NorthEast');
legend([p1 p2 p3 p4],'True Values', 'Measurements', 'Estimate','95% Confidence Interval')
%title('Building Height - 95% Confidence Interval of the Estimates', 'FontSize', 12, 'FontWeight', 'Bold');
