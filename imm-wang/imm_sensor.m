function [sensor_mu, sensor_var] = imm_sensor(imm_labels, sensor_data)

  % Get the dimensions (which is the max in the labels)
%  m = max(imm_labels);
  m = 2;

  % Find the sensor data mean according to the labels
  sensor_mu = nan(m,1);
  sensor_var = nan(m,1);
  for k = 1:m
    sensor_mu(k) = mean(sensor_data(find(imm_labels == k)));
    sensor_var(k) = var(sensor_data(find(imm_labels == k)));
  end

end %EOF
