function plot_velocity_triangles_angles(v1, v2, u1, u2, w1, w2, chord_length, inclination_angle, beta1, beta2, alpha1, alpha2)
    % Plot w1, u1 with the tail of u1 on the head of w1, and v1 connecting tail of w1 to head of u1
    figure;
    hold on;
    axis equal;
    grid on;

    v1 = [v1 * cos(alpha1), v1 * sin(alpha1)];
    v2 = [v2 * cos(alpha2), v2 * sin(alpha2)];
    u1 = [0, u1];
    u2 = [0, u2];
    w1 = [w1 * cos(beta1), w1 * sin(beta1)];
    w2 = [w2 * cos(beta2), w2 * sin(beta2)];

    mag_v1 = norm(v1);
    mag_v2 = norm(v2);
    mag_u1 = norm(u1);
    mag_u2 = norm(u2);
    mag_w1 = norm(w1);
    mag_w2 = norm(w2);
    
    % i rescale here so that the diagram is readable
    max_magnitude = max([mag_v1, mag_v2, mag_u1, mag_u2, mag_w1, mag_w2]);
    scale_factor = 2 * chord_length / max_magnitude;
    
    % Scale all vectors
    v1 = v1 * scale_factor;
    v2 = v2 * scale_factor;
    u1 = u1 * scale_factor;
    u2 = u2 * scale_factor;
    w1 = w1 * scale_factor;
    w2 = w2 * scale_factor;
    
    quiver(0, 0, w1(1), w1(2), 0, 'r', 'LineWidth', 2);
    quiver(w1(1), w1(2), u1(1), u1(2), 0, 'b', 'LineWidth', 2);
    quiver(0, 0, v1(1), v1(2), 0, 'g', 'LineWidth', 2);
    
    theta = inclination_angle;
    chord_vec = chord_length * [cos(theta), sin(theta)];
    quiver(w1(1), w1(2), chord_vec(1), chord_vec(2), 0, 'k', 'LineWidth', 2);
    
    chord_end = w1 + chord_vec;
    
    quiver(chord_end(1), chord_end(2), w2(1), w2(2), 0, 'r', 'LineWidth', 2);
    quiver(chord_end(1) + w2(1), chord_end(2) + w2(2), u2(1), u2(2), 0, 'b', 'LineWidth', 2);
    quiver(chord_end(1), chord_end(2), v2(1), v2(2), 0, 'g', 'LineWidth', 2);
    
    text(w1(1)/2, w1(2)/2, 'w1', 'FontSize', 12);
    text(w1(1) + u1(1)/2, w1(2) + u1(2)/2, 'u1', 'FontSize', 12);
    text(v1(1)/2, v1(2)/2, 'v1', 'FontSize', 12);
    text(chord_end(1) + w2(1)/2, chord_end(2) + w2(2)/2, 'w2', 'FontSize', 12);
    text(chord_end(1) + w2(1) + u2(1)/2, chord_end(2) + w2(2) + u2(2)/2, 'u2', 'FontSize', 12);
    text(chord_end(1) + v2(1)/2, chord_end(2) + v2(2)/2, 'v2', 'FontSize', 12);
    
    xlabel('X');
    ylabel('Y');
    title('Velocity Triangles with Chord Inclination');
    hold off;
end
