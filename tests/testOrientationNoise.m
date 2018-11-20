figure, hold on
rng(1)
s = 0;
n = 0;
dT = 0.035/0.33/16;
omega = 2*pi*0.6;
amplitute = pi/4;
t = 0;
T = 5;
while t<T
s = s + amplitute*omega*cos(omega*t)*dT;
n = n + amplitute*omega*cos(omega*t)*dT + 2*randn()*(dT);
t = t + dT;
plot(t,s,'rx')
plot(t,n,'bx')
end