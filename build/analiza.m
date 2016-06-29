r = rs;
v = vs;
a = as;

t = rs(:,1);
plot( r(:,2), 'b*');
hold on;
plot( v(:,1), 'g*');
plot( a(:,1), 'r*');
%legend( 'position', 'velocity', 'acceleration' )

plot( sin(t), 'k')
plot( cos(t), 'g')
plot( -sin(t),'r')

hold off;
