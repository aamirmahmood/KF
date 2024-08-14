clc, close all, clear all;

% Measurements vector
z = [30221 30453 30906 30999 31368 31978 32526 33379 34698 36275];

% Filter Parameters
alpha = 0.5; 
beta  = 0.4;
gamma = 0.1;

% Tracking 
delta_t = 5; % Track cycle

% Initialization 
range_estimate        = 30000; % in meters 
velocity_estimate     = 50;
acceleration_estimate = 0;

% Prediction at Iteration Zero
range_estimate_pred        = zeros(1,length(z));
velocity_estimate_pred     = zeros(1,length(z));
acceleration_estimate_pred = zeros(1,length(z));

range_estimate_pred(1,1)        = range_estimate + (velocity_estimate * delta_t) + (0.5*acceleration_estimate*delta_t^2);
velocity_estimate_pred(1,1)     = velocity_estimate + (acceleration_estimate*delta_t);
acceleration_estimate_pred(1,1) = acceleration_estimate; 

for ii = 1:length(z)
    
    % State Updates
    range_estimate_current(ii)        = range_estimate_pred(1,ii) + alpha*(z(ii) - range_estimate_pred(1,ii));
    velocity_estimate_current(ii)     = velocity_estimate_pred(1,ii) + beta*((z(ii) - range_estimate_pred(1,ii))/delta_t);
    acceleration_estimate_current(ii) = acceleration_estimate_pred(1,ii) + gamma*((z(ii) - range_estimate_pred(1,ii))/(0.5*delta_t^2));
    
    % Prediction Updates
    range_estimate_pred(1, ii+1)       = range_estimate_current(ii) + (delta_t * velocity_estimate_current(ii)) + (0.5*acceleration_estimate_current(ii)*delta_t^2);
    velocity_estimate_pred(1,ii+1)     = velocity_estimate_current(ii) + (delta_t*acceleration_estimate_current(ii));
    acceleration_estimate_pred(1,ii+1) = acceleration_estimate_current(ii);
end

% Time
time = 0:5:(length(z)-1)*5;

% True Value
% true_range = range_estimate + ((time+5) * velocity_estimate); 




%% Baseline
s_init = 30200;
t_15   = 0:1:15;
t_35   = 16:1:50;
v      = 50; 
a      = 8;

delta_t = t_15(2) - t_15(1);

% s_15 = zeros(1,length()

for ii=1:length(t_15)
    s_15(ii)  = s_init + (v*delta_t);
    s_init    = s_15(ii);
end

s_init = s_15(end); 
v_init = v;
for ii=1:length(t_35)
    s_35(ii) = s_init + (v_init*delta_t) + (0.5*a*delta_t^2);
    v_35(ii) = v_init + (a*delta_t);
    v_init = v_35(ii);
    s_init = s_35(ii);
end
t = [t_15 t_35];
s = [s_15 s_35];

%%

figure;
subplot(3,1,1);
% plot(time, true_range, '-*', 'LineWidth', 2);
plot(t(1:end-5), s(1:end-5), '-*', 'LineWidth', 2);
hold on;
plot(time, z, '-s', 'LineWidth', 2);
plot(time, range_estimate_current, '-o', 'LineWidth', 2);
plot(time, range_estimate_pred(1:end-1), '-d', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
title('Range vs. Time');
legend('True Values', 'Measurements', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');

% figure;
subplot(3,1,2);

plot(t(1:end-5), [ones(1,(length(t_15)))*v v_35(1:end-5)], '-*', 'LineWidth', 2);
hold on;
plot(time, velocity_estimate_current, '-o', 'LineWidth', 2);
plot(time, velocity_estimate_pred(1:end-1), '-o', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Velocity vs. Time')
grid on;
legend('True Values', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');

% figure;
subplot(3,1,3);
plot(t, [zeros(1,(length(t_15))) ones(1,length(t_35))*a], '-k', 'LineWidth', 2);
hold on;
plot(time, acceleration_estimate_current, '-o', 'LineWidth', 2);
plot(time, acceleration_estimate_pred(1:end-1), '-o', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Acceleration (m/s^2)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Acceleration vs. Time')
grid on;
legend('True Values', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');

