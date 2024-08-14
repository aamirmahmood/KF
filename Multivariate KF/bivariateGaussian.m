clc;clear all;close all;
% Define the grid of values for X and Y
[x, y] = meshgrid(-3:.1:3, -3:.1:3);

% Define the parameters for the means and standard deviations
mu_x = 0;
mu_y = 0;
sigma_x = 1;
sigma_y = 1;

% Define the correlation coefficient
rho = 0.8; % Correlation coefficient

% Compute the covariance matrix elements
sigma_xy = rho * sigma_x * sigma_y;
covariance_matrix = [sigma_x^2, sigma_xy; sigma_xy, sigma_y^2];

% Compute the inverse and determinant of the covariance matrix
inv_covariance_matrix = inv(covariance_matrix);
det_covariance_matrix = det(covariance_matrix);

% Calculate the Z values for the correlated bivariate Gaussian distribution
z = zeros(size(x));
for i = 1:numel(x)
    diff_vec = [x(i); y(i)] - [mu_x; mu_y];
    z(i) = (1 / sqrt((2*pi)*det_covariance_matrix)) * exp(-0.5 * (diff_vec' * inv_covariance_matrix * diff_vec));
end

% Reshape Z to the size of X and Y
z = reshape(z, size(x));

% Create the 3D surface plot
surf(x, y, z)
colormap('hot')
%shading interp
xlabel('X', 'FontSize',12,'FontWeight','bold')
ylabel('Y', 'FontSize',12,'FontWeight','bold')
zlabel('Z', 'FontSize',12,'FontWeight','bold')

% Adjust the axis limits and labels
xlim([-3 3])
ylim([-3 3])
zlim([0 max(z, [], 'all')]) % Use the maximum value of Z for the upper Z limit


% Add contour lines to the plot
hold on
contour3(x, y, z); % Draw 10 contour lines in black

% Set the view angle
view(-20, 20)

% Save the figure to a file
saveas(gcf, 'correlated_gaussian_surface.eps')



% Define the mean and covariance matrix
mu = [0 1];
rho = 0.5;
Sigma = [1 rho; rho 1];

% Create a grid of (x,y) points
[x, y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

% Calculate the bivariate Gaussian values
z = mvnpdf([x(:) y(:)], mu, Sigma);
z = reshape(z, size(x));

% Plot the surface
figure;
surf(x, y, z);

% Customize the plot
xlabel('X', 'FontSize', 14, 'Color', 'red');
ylabel('Y', 'FontSize', 14, 'Color', 'red');
zlabel('Z', 'FontSize', 14, 'Color', 'red');
colorbar;
title('Bivariate Gaussian Distribution');

% Adjust the view angle for better visualization
view(60, 30);

