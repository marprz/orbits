a = sat1d1; 
b = sat1d2;
c = sat1d3;

d = [ a; b; c ];
tmax = 143.5;
t = 0:0.5:tmax;
subplot(3,1,1)
plot(t,a(:,1)-b(:,1),'g');
hold on;
plot(t,a(:,1)-c(:,1),'b');

n = 288;
e = a(5:288,1)-b(1:284,1)
