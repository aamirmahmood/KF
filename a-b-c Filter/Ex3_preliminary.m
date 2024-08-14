% Example 3

clc;close all; clear all;
v = 40; 

s_init = 30200;

t       = 0:0.1:50;
delta_t = t(2)-t(1);

for ii=1:length(t)
    s_new(ii)  = s_init + v*delta_t;
    s_init = s_new(ii);
end


subplot(2,1,1)
plot(t, s_new*1e-3, '-r', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Range vs. Time')
grid on;

subplot(2,1,2)
plot(t, v*ones(1,length(t)), '-b', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Velocity vs. Time')
grid on;


%%
% Assume, a fighter aircraft moving at a constant velocity of 50\,m/s for 15 seconds. 
% Then the aircraft accelerates with a constant acceleration of $8\,m/s^2$ for another 35 seconds.
clear all;

s_init = 30200;
t_15   = 0:0.1:15;
t_35   = 15.1:0.1:70;
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

% t = [t_15 t_35];

% s_init = s_15(end); 
% s_init = 30200;
% for ii=1:length(t)
%     s_35(ii) = s_init + (v*t) + (0.5 * a * t^2);
%     s_init = s_35(ii);
% end

% s_15 = s_init + (v*t_15);
% s_35 = 30200 + (v*t) + (0.5*a*t.^2);

figure;

t = [t_15 t_35];
s = [s_15 s_35];

subplot(3,1,1)
plot(t, s*1e-3, '-r', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Range (m)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Range vs. Time')
grid on;

subplot(3,1,2)
plot(t, [ones(1,(length(t_15)))*v v_35], '-b', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Velocity vs. Time')
grid on;


subplot(3,1,3)
plot(t, [zeros(1,(length(t_15))) ones(1,length(t_35))*a], '-k', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'Bold');
title('Acceleration vs. Time')
grid on;





