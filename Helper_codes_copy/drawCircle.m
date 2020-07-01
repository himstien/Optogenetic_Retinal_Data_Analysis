function drawCircle(center, radius, color)

x = center(1);
y = center(2);

ang = 0:0.01:2*pi;

xp = radius*cos(ang);
yp = radius*sin(ang);

plot(x+xp, y+yp, color, 'LineWidth', 2);