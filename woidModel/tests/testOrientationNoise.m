figure, hold on
s = 0;
n = 0;
for t=1:1000
s = s + 0.1*cos(0.1*t);
n = n + 0.1*cos(0.1*t) + 0.1*randn();
plot(t,s,'rx')
plot(t,n,'bx')
end