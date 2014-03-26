[x,y,z] = meshgrid(1:0.1:2, -2:.2:2, 0:0.1:1);
C = 1/x
figure
surf(x,y,z,C);
colormap(cool)
axis equal