clc, close all, clear all;


% Measurements vector
z = [30110 30265 30740 30750 31135 31015 31180 31610 31960 31865];

% Filter Parameters
alpha = 0.8; 
beta  = 0.5;

% Tracking 
delta_t = 5; % Track cycle

% Initialization 
range_estimate    = 30000; % in meters 
velocity_estimate = 40;

% Prediction at Iteration Zero
range_estimate_pred = zeros(1,length(z));
velocity_estimate_pred = zeros(1,length(z));
range_estimate_pred(1,1)    = range_estimate + (delta_t * velocity_estimate);
velocity_estimate_pred(1,1) =  velocity_estimate;


for ii = 1:length(z)
    
    % State Updates
    range_estimate_current(ii)    = range_estimate_pred(1,ii) + alpha*(z(ii) - range_estimate_pred(1,ii));
    velocity_estimate_current(ii) = velocity_estimate_pred(1,ii) + beta*((z(ii) - range_estimate_pred(1,ii))/delta_t);
    
    % Prediction Updates
    range_estimate_pred(1, ii+1) = range_estimate_current(ii) + (delta_t * velocity_estimate_current(ii));
    velocity_estimate_pred(1,ii+1) = velocity_estimate_current(ii);
end

% Time
time = 0:5:(length(z)-1)*5;

% True Value
true_range = range_estimate + ((time+5) * velocity_estimate); 


plot(time, true_range, '-*', 'LineWidth', 2);
hold on;
plot(time, z, '-s', 'LineWidth', 2);
plot(time, range_estimate_current, '-o', 'LineWidth', 2);
plot(time, range_estimate_pred(1:end-1), '-d', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('True Values', 'Measurements', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');



