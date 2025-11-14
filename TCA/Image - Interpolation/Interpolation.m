clear;clc;
filename = 'up_down.mat'; % Set the contour filename from grabit
output_filename = 'up_down.csv'; % Change for your desired output filename
n_points = 100; % Set the number of points you wat to have.
axis_array = [0,0.6,0,0.2]; %Set your axis for the plot


%%%%%%%%%%%%%%% SCRIPT %%%%%%%%%%%%%%%
S  = load(filename); fn = fieldnames(S);
M  = S.(fn{1});
x = M(:,1);
y = M(:,2);

% build interpolant
f = make_interpolant(x, y);
xq = linspace(min(x), max(x), n_points);
yq = f(xq);

data = [xq(:), yq(:)];      % make sure they are columns
writematrix(data, output_filename);

% --- plot ---
figure;
plot(x, y, 'o', 'DisplayName', filename); hold on;
plot(xq, yq, '-', 'DisplayName', 'interpolant');
xlabel('x'); ylabel('y');
legend('Location','best','Interpreter','none');
grid on;axis equal;   
axis(axis_array);
