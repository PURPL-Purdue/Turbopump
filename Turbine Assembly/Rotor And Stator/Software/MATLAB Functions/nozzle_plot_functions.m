function funs = nozzle_plot_functions
    % function that returns a struct of function handles
    funs.plot_nozzle = @plot_nozzle;
    funs.exp_scale = @exp_scale;
    funs.save_to_stl = @save_to_stl;
end


% function definitions
function [X,Y,Z] = plot_nozzle(A_inlet, A_throat, A_exit, inlet_len, outlet_len)
    figure;
    num_points = 100;
    r_inlet = sqrt(A_inlet/pi);
    r_throat = sqrt(A_throat/pi);
    r_exit = sqrt(A_exit/pi);
    % outlet_len = (r_exit - r_throat) * tand(angle);
    throat_len = 0.1 * inlet_len;
    x_total = inlet_len + throat_len + outlet_len;
    x_step = x_total / num_points;
    [X, Thetas] = meshgrid(0:x_step:x_total, 0:(2*pi/1000):2*pi);
    R = (r_inlet + (r_throat - r_inlet) .* exp_scale(X ./ inlet_len)) .* (X <= inlet_len);
    R = R + r_throat .* (inlet_len < X & X <= (inlet_len + throat_len));
    R = R + (r_throat + ((r_exit - r_throat)/outlet_len .* (X - (inlet_len + throat_len)) )) .* ((inlet_len + throat_len) < X & X <= x_total);
    Y = R .* cos(Thetas);
    Z = R .* sin(Thetas);
    C = X;

    surf(X,Y,Z,C)
    shading interp
    axis equal
    xlabel("length [m]")
    ylabel("width [m]")
    zlabel("height [m]")
    title("Plot of nozzle geometry")
    colorbar

    minZ = min(Z(:));
    maxZ = max(Z(:));  

    figure;
    contour(X,Y,Z,linspace(minZ, maxZ / 2, 10))
    xlabel("length [m]")
    ylabel("width [m]")
    title("Contour Plot for nozzle")
end

function percentage=exp_scale(x)
    percentage = (exp(x) - 1) ./ (exp(1) - 1);
end

function save_to_stl(X, Y, Z, filename)
    [faces, vertices] = surf2patch(X, Y, Z, 'triangles');
    tr = triangulation(faces, vertices);
    stlwrite(tr, filename);    
    disp(['Surface exported to ', filename]);
end
