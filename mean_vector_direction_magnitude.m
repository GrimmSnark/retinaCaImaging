function stats = mean_vector_direction_magnitude(vector_input)
% Average vectors and output both mean angle and magnitude.
% input is rows of [angle_in_degrees magnitude_at_that_angle]

x_values = []; y_values = []; % values of vector components in the two dimensions 
for i=1:size(vector_input,1)
    x_values = [x_values; vector_input(i, 2)*cosd(vector_input(i, 1))];
    y_values = [y_values; vector_input(i, 2)*sind(vector_input(i, 1))];
    
end

% calculate vector from means
mean_angle_raw = atan2(mean(y_values), mean(x_values)); % this is raw degree value based on atan
mean_angle_raw = mean_angle_raw .* (mean_angle_raw >= 0) + (mean_angle_raw + 2 * pi) .* (mean_angle_raw < 0); % correct to give full range from 0 -> 2pi
stats.mean_angle_degrees = rad2deg(mean_angle_raw); % convert to degrees

if mean(abs(y_values)) > mean(abs(x_values)) % try to avoid dividing by zero and getting NaN
    stats.mean_magnitude = abs(mean(y_values)/sind(stats.mean_angle_degrees));
else
    stats.mean_magnitude = abs(mean(x_values)/cosd(stats.mean_angle_degrees));
end

end

