T = readtable("nozzle_contour.csv");
x = T.('x');
y = T.('y');

plot(x,y);
axis equal

A = pi*(y/1000).^2;