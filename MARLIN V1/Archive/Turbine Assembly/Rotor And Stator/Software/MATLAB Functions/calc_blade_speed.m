function U = calc_blade_speed(radius, rpm)
    % rpm to angular velocity (rad/s)
    omega = rpm * 2 * pi / 60;
    U = [radius * omega, 0];  % blade speed (m/s)
end
