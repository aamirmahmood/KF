clc, close all, clear all;

h_true = 50;
% Measurements vector
z_1 = [49.03, 48.44, 55.21, 49.98, 50.6, 52.61, 45.87, 42.64, 48.26, 55.84];
z_2 = [49.3 49.1 45.2 45.4 49.7 54.6 48.03 55.85 51.95 55.15 52.01 43.20 50.7 47.7 50.96 54.8 41 50 52.98 ...
       54.2 54.85 49.2 48 45.5 44.98 48 61.2 48.15 56.15 57.2 42.05 48 64.2 50 59.91 45.65 41.25 62.6 51.5 46.5];

z = [z_1 z_2];


% Measurement uncertainty 
sigma_std = 2; 
r = sigma_std^2;


%% ITERATION ZERO

% Initialization 

h_estimate = 60;                % Estimated build height
uncertainity_estimate = 225;    % Estimate uncertainty with std. σ = 15m

% Prediction at Iteration Zero
h_estimate_pred(1,1) = h_estimate;   % Since Dynamic Model is constant
uncertainity_estimate_pred(1,1) = uncertainity_estimate; % The extrapolated estimate uncertainty (variance) also doesn’t change

for ii = 1:length(z)
    
%   1) Measure
    z_n = z(ii); % New Measurement
    r_n = r;    % Measurement uncertainty
    
    
    % Kalman gain
    K(ii) = uncertainity_estimate_pred(1,ii) / (uncertainity_estimate_pred(1,ii) + r_n);
    
    
%   3) State Update
    h_estimate_current(ii)    = h_estimate_pred(1,ii) + K(ii)*(z(ii) - h_estimate_pred(1,ii)); %Estimating the current state
    uncertainity_estimate_current(ii) = (1-K(ii))* uncertainity_estimate_pred(1,ii); %Current estimate uncertainty
    
    % Prediction Updates
    h_estimate_pred(1, ii+1) = h_estimate_current (ii);
    uncertainity_estimate_pred(1,ii+1) = uncertainity_estimate_current(ii);
end

x = 1:length(z);

figure();

plot(x, ones(1,length(z))*h_true, '-d', 'LineWidth', 2);
hold on;
plot(x, z, '-*', 'LineWidth', 2);

plot(x, h_estimate_current, '-s', 'LineWidth', 2);
plot([0], h_estimate, 'd', 'LineWidth', 4);

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Height (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('True Values', 'Measurements', 'Estimates', 'Init', 'FontSize', 12, 'Location', 'NorthEast');
title('Building Height', 'FontSize', 12, 'FontWeight', 'Bold');


%%
figure;

p1 = plot(x, ones(1,length(z))*h_true, '-d', 'LineWidth', 2);
hold on;
p2= plot(x, h_estimate_current, '-s', 'LineWidth', 2);
% 95% Confidence Interval

CI = 95;
temp = (1 - CI/100)/2;

low = h_estimate_current - norminv(temp)*sqrt(uncertainity_estimate_current);
up  = h_estimate_current + norminv(temp)*sqrt(uncertainity_estimate_current);

plot(x, low, 'r', 'LineWidth',1);
plot(x, up, 'r', 'LineWidth',1);
x2 = [x, fliplr(x)];
inBetween = [low, fliplr(up)];
p3 = patch(x2, inBetween, 'y');
%legend(p1,'Area 1')
alpha(0.25)
grid on;

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Height (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
%legend('True Values', 'Estimates', 'FontSize', 12, 'Location', 'NorthEast');
legend([p1 p2 p3],'True Values','Estimates','95% Confidence Interval')
title('Building Height - 95% Confidence Interval of the Estimates', 'FontSize', 12, 'FontWeight', 'Bold');

%%
figure;

plot(x, ones(1,length(z))*r, '-d', 'LineWidth', 2);
hold on;
plot(x, uncertainity_estimate_current, '-*', 'LineWidth', 2);

% plot(x, h_estimate_current, '-s', 'LineWidth', 2);
plot([0], uncertainity_estimate, 'd', 'LineWidth', 4);

xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Uncertainty (m^2)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('Measurement Uncertainty', 'Estimates Uncertainty', 'Init', 'FontSize', 12, 'Location', 'NorthEast');
title('Uncertainties', 'FontSize', 12, 'FontWeight', 'Bold');

figure;

plot(x, K, '-d', 'LineWidth', 2);
xlabel('Measurement Number', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
title('Kalman Gain', 'FontSize', 12, 'FontWeight', 'Bold');


