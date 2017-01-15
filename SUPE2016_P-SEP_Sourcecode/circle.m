function circle(x,y,r)
t = (0:0.01:2*pi)';
xp = r*cos(t);
yp = r*sin(t);
xr=x+xp;
yr=y+yp;
plot(xr,yr);
end
