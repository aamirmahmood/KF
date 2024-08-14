clc; close all; clear all;

d = 300; 
measurement_angle = 0:1:90; % measurement angle in degrees
altitude = d*tand(measurement_angle);

plot (measurement_angle, altitude, 'b-', 'LineWidth', 2)


% select a time index

angle_selected = 83;    
altitude_selected = d*tand(angle_selected);

angle_uncertainty_std = 1;
angle_uncertainty = angle_selected+angle_uncertainty_std*randn(1,1000000);
altitude_uncertainty = d*tand(angle_uncertainty);

subplot(2,2,2)
plot (measurement_angle,altitude, 'b-', 'LineWidth', 2)
hold on;
line([angle_selected angle_selected], [altitude(1) altitude_selected], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
line([measurement_angle(1) angle_selected], [altitude_selected altitude_selected], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel('Angle Measurement (deg)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('Altitude (m)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
grid on;


subplot(2,2,4)
%histogram(t_uncertainty,'Normalization','pdf')
[f, xi] = ksdensity(angle_uncertainty(1:end));
plot(xi, f, 'b', 'LineWidth', 2);
xlabel('Angle uncertainty (deg)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
ylabel('PDF', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0, 0]);
line([angle_selected angle_selected], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
xlim([0, max(xi)]); % Set the x-axis limits to start from 0

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
