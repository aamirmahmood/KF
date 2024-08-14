%In city A, the mean delivery time is 30 minutes, and the standard deviation is 5 minutes.

mu    = 30; % mean delivery time in minutes
sigma = 5;  % standard deviation in minutes

%% Q1: What is the likelihood that the pizza in city A would be delivered within 33 minutes?

% Answer: find p(0 <= x <= 33)


x_up  = 33;
x_low = 0;

% Method 1: Calculate the z score

z_up  = (x_up - mu)/sigma;
z_low = (x_low - mu)/sigma;

% Likelihood of receiving pizaa within x_up minsutes
LH = normcdf(z_up) - normcdf(z_low)

% Method 2
LH = normcdf(33, 30, 5)


%% What is the 80th percentile for the pizza delivery time in the city A?

% Method 1: find the z-score using the cumulative distribution value closest to 0.8.
% x = z*sigma + mu

% Method 2: 

DeliveryTime_80Percentil = norminv(0.8, 30, 5)

%% Confidence Interval

% assuming 90% confidence interval
CI = 90;
temp = (1 - CI/100)/2;

range_l = mu + norminv(temp)*sigma
range_u = mu - norminv(temp)*sigma
