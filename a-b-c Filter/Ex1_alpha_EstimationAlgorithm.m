clc, close all, clear all;


% Measurements vector
z = [1030 989 1017 1009 1013 979 1008 1042 1012 1011];

% Initialization 
x_estimate      = 1000; 
x_estimate_pred = x_estimate;

for ii = 1:length(z)
    alpha(ii) = 1/ii; % Kalman gain
    x_estimate_current(ii) = x_estimate_pred + alpha(ii) * (z(ii) - x_estimate_pred); % state update equation
    x_estimate_pred = x_estimate_current(ii); % the current estimate becomes the prior estimate in the next iteration
end

% True Value
x = ones(1,length(z))*1010;

plot(1:length(z), x, '-*', 'LineWidth', 2);
hold on;
plot(1:length(z), z, '-s', 'LineWidth', 2);
plot(1:length(z), x_estimate_current, '-o', 'LineWidth', 2);


xlabel('Iterations', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Weight (g)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('True Values', 'Measurements', 'Estimates', 'FontSize', 12, 'Location', 'NorthWest');



