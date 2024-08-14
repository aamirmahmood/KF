function r_ellipses = calculate_error_ellipses(x_y_estimate, x_y_uncertainty, CI)

    % Chi-square value for the 95% confidence interval
    chisquare_val = sqrt(chi2inv(CI/100, 2));

    %Initialize cell array to store ellipses
    r_ellipses = cell(1, length(x_y_estimate));

    for ii = 1:length(x_y_estimate)
        % Calculate the eigenvectors and eigenvalues
        covariance = x_y_uncertainty(:, :, ii);
        [eigenvec, eigenval] = eig(covariance);

        % Get the index of the largest eigenvector
        [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
        largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

        % Get the largest eigenvalue
        largest_eigenval = max(max(eigenval));

        % Get the smallest eigenvector and eigenvalue
        if largest_eigenvec_ind_c == 1
            smallest_eigenval = max(eigenval(:, 2));
            smallest_eigenvec = eigenvec(:, 2);
        else
            smallest_eigenval = max(eigenval(:, 1));
            smallest_eigenvec = eigenvec(:, 1);
        end

        % Calculate the angle between the x-axis and the largest eigenvector
        angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

        % Shift angle to be between 0 and 2pi
        if angle < 0
            angle = angle + 2 * pi;
        end

        % Get the coordinates of the data mean
        X0 = x_y_estimate(1, ii);
        Y0 = x_y_estimate(2, ii);

        % Get the 95% confidence interval error ellipse
        a = chisquare_val * sqrt(largest_eigenval);
        b = chisquare_val * sqrt(smallest_eigenval);

        % Define the ellipse in x and y coordinates
        theta_grid = linspace(0, 2 * pi);
        ellipse_x_r = a * cos(theta_grid);
        ellipse_y_r = b * sin(theta_grid);

        % Define a rotation matrix
        R = [cos(angle) sin(angle); -sin(angle) cos(angle)];

        % Rotate the ellipse to the angle phi
        r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;
        
        r_ellipse = [r_ellipse(:, 1) + X0, r_ellipse(:, 2) + Y0];
        % Draw the error ellipse
        %plot(r_ellipse(:, 1) + X0, r_ellipse(:, 2) + Y0, '-');

        % Store the ellipse
        r_ellipses{ii} = r_ellipse;
    end

%     hold off;
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     title('95% Confidence Interval Error Ellipses');
%     axis equal;
%     grid on;
end
