close all
% Parameters for the bivariate normal distribution
mu = [0 0]; % Mean
sigma = [1 0.8; 0.8 1]; % Covariance matrix

% Create grid and multivariate normal values
[x, y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));
X = [x(:) y(:)];
Z = mvnpdf(X, mu, sigma);
Z = reshape(Z, size(x));


% Plot the surface
figure;
surf(x, y, Z, 'FaceAlpha', 0.8);
colormap parula
hold on
alpha 0.5;

% Add contours at specific Z-level (e.g., Z = -0.1 for visibility)
contourf(x, y, Z, 15,'LineWidth', 2); % Draw 20 contour lines in black
zlim([min(Z(:)) max(Z(:))]);
hold off

xlabel('X');
ylabel('Y');
zlabel('Probability Density (Z)');
%title('Bivariate Normal Distribution with Contours');

view(-20, 20)

% Save the figure to a file
%saveas(gcf, 'correlated_gaussian_surface_rh0.8.eps')
