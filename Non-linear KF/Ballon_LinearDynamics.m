clc; close all; clear all;

%% Method 1: Constant acceleration

% Given parameters
initial_altitude = 308; % meters
final_altitude = 315; % meters
initial_time = 2.0e-6; % seconds
final_time = 2.1e-6; % seconds
dt = 1e-9; % Time step

% Compute the time interval
time_interval = final_time - initial_time;

% Compute the required constant acceleration
acceleration = 2 * (final_altitude - initial_altitude) / (time_interval^2);

% Time vector for simulation
t = initial_time:dt:final_time;
num_steps = length(t);

% Initialize altitude and velocity
altitude = zeros(1, num_steps);
velocity = zeros(1, num_steps);
altitude(1) = initial_altitude;
velocity(1) = 0;

% Simulate the balloon motion
for k = 2:num_steps
    velocity(k) = velocity(k-1) + acceleration * dt;
    altitude(k) = altitude(k-1) + velocity(k-1) * dt + 0.5 * acceleration * dt^2;
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(t, altitude, 'b', 'DisplayName', 'Altitude');
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Balloon Altitude Simulation');
legend;

subplot(2, 1, 2);
plot(t, velocity, 'r', 'DisplayName', 'Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Balloon Velocity Simulation');
legend;


%% Method 2: Constant velocity
% Given parameters
initial_altitude = 308; % meters
final_altitude = 315; % meters
initial_time = 2.0e-6; % seconds
final_time = 2.1e-6; % seconds
dt = 1e-9; % Time step

% Compute the time interval
time_interval = final_time - initial_time;

% Compute the required constant velocity
velocity = (final_altitude - initial_altitude) / time_interval;

% Time vector for simulation
t = initial_time:dt:final_time;
num_steps = length(t);

% Initialize altitude
altitude = zeros(1, num_steps);
altitude(1) = initial_altitude;

% Simulate the balloon motion with constant velocity
for k = 2:num_steps
    altitude(k) = altitude(k-1) + velocity * dt;
end

% Plot the results
figure;
plot(t, altitude, 'b', 'DisplayName', 'Altitude');
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Balloon Altitude Simulation with Linear Increase');
legend;

%% Let's do our own as many things are unknow to create the figure like in the book
clear all;
c = 3.8e9;

t = 2.0e-6:0.001e-6:2.1e-6;
altitude = (c/2) * t;

% select a time index

t_selected = t(find(t==2.05e-6));
altitude_selected = (c/2) * t_selected;

t_uncertainty_std = 0.01e-6;
t_uncertainty = t_selected+t_uncertainty_std*randn(1,1000000);
altitude_uncertainty = (c/2) * t_uncertainty;

subplot(2,2,2)
plot (t,altitude, 'b-', 'LineWidth', 2)
hold on;
line([t_selected t_selected], [altitude(1) altitude_selected], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
line([t(1) t_selected], [altitude_selected altitude_selected], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel('Time measurement (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Altitude (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;


subplot(2,2,4)
%histogram(t_uncertainty,'Normalization','pdf')
[f, xi] = ksdensity(t_uncertainty);
plot(xi, f, 'b', 'LineWidth', 2);
xlabel('Time uncertainty (s)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('PDF', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
line([t_selected t_selected], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

grid on;

subplot(2,2,1)
[f, xi] = ksdensity(altitude_uncertainty);
plot(xi,f, 'r', 'LineWidth', 2);
xlabel('PDF', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Altitude uncertainty (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;
hold on;
line([altitude_selected altitude_selected], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

set(gca, 'View', [-90 90])
%line(xlim, [altitude_selected altitude_selected], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
