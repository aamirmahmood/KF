clc, close all, clear all;

T_true = 50;
z_true = [50.479 51.025 51.5 52.003 52.494 53.002 53.499 54.006 54.498 54.991];


% Measurements vector
z = [50.45 50.967 51.6 52.106 52.492 52.819 53.433 54.007 54.523 54.99];

% Measurement uncertainty 
sigma_std = 0.1; 
r = sigma_std^2;

% process noise variance
q = 0.15;

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

subplot(2,1,1);
p1 = plot(x, z_true, '-d', 'LineWidth', 2);
hold on;
p2 = plot(x, z, '-*', 'LineWidth', 2);

p3 = plot(x, T_estimate_current, '-s', 'LineWidth', 2);
% plot([0], T_estimate, 'd', 'LineWidth', 4);

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
ylabel('Temperature (C)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
%legend('True Values', 'Measurements', 'Estimate', 'FontSize', 12, 'Location', 'NorthEast');
legend([p1 p2 p3 p4],'True Values', 'Measurements', 'Estimate','95% Confidence Interval')
title('Temperature', 'FontSize', 12, 'FontWeight', 'Bold');

% 
% subplot(3,1,2);
% % plot(x, ones(1,length(z))*r, '-d', 'LineWidth', 2);
% % hold on;
% plot(x, uncertainity_estimate_current, '-*', 'LineWidth', 2);
% 
% % plot(x, h_estimate_current, '-s', 'LineWidth', 2);
% % plot([0], uncertainity_estimate, 'd', 'LineWidth', 4);
% 
% xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('Uncertainty', 'FontSize', 12, 'FontWeight', 'Bold');
% grid on;
% legend('Estimate Uncertainty', 'FontSize', 12, 'Location', 'NorthEast');
% title('Uncertainty', 'FontSize', 12, 'FontWeight', 'Bold');

subplot(2,1,2);

plot(x, K, '-d', 'LineWidth', 2);
xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
title('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');


