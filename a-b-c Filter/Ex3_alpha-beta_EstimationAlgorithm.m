    clc, close all, clear all;

% Measurements vector
z = [30160 30365 30890 31050 31785 32215 33130 34510 36010 37265];
z = [30221 30453 30906 30999 31368 31978 32526 33379 34698 36275];

% Filter Parameters
alpha = 0.2; 
beta  = 0.1;

% Tracking 
delta_t = 5; % Track cycle in seconds

% Initialization 
range_estimate    = 30000; % in meters 
velocity_estimate = 50;

% Prediction at Iteration Zero
range_estimate_pred = zeros(1,length(z));
velocity_estimate_pred = zeros(1,length(z));
range_estimate_pred(1,1) = range_estimate + (delta_t * velocity_estimate);
velocity_estimate_pred(1,1)   =  velocity_estimate;


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
% true_range = range_estimate + ((time+5) * velocity_estimate); 




%% Baseline

s_init = 30200;
t_15   = 0:1:15;
t_35   = 16:1:50;
v      = 50; 
a      = 8;

delta_t = t_15(2) - t_15(1);


for ii=1:length(t_15)
    s_15(ii)  = s_init + (v*delta_t);
    s_init = s_15(ii);
end

s_init = s_15(end); 
v_init = v;
for ii=1:length(t_35)
    s_35(ii) = s_init + (v_init*delta_t) + (0.5*a*delta_t^2);
    v_35(ii) = v_init + (a*delta_t);
    v_init = v_35(ii);
    s_init = s_35(ii);
end

% figure;

t = [t_15 t_35];
s = [s_15 s_35];

% plot(t, s, '-g', 'LineWidth', 2);
% xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
% title('Range vs. Time')
% grid on;

%%

% plot(time, true_range, '-*', 'LineWidth', 2);
plot(t(1:end-5), s(1:end-5), '-*', 'LineWidth', 2);
hold on;
plot(time, z, '-s', 'LineWidth', 2);
plot(time, range_estimate_current, '-o', 'LineWidth', 2);
plot(time, range_estimate_pred(1:end-1), '-d', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
grid on;
legend('True Values', 'Measurements', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');

figure;

plot(t(1:end-5), [ones(1,(length(t_15)))*v v_35(1:end-5)], '-*', 'LineWidth', 2);
hold on;
plot(time, velocity_estimate_current, '-o', 'LineWidth', 2);
plot(time, velocity_estimate_pred(1:end-1), '-o', 'LineWidth', 2);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'Bold');
% title('Velocity vs. Time')
grid on;
legend('True Values', 'Estimates', 'Prediction', 'FontSize', 12, 'Location', 'NorthWest');

