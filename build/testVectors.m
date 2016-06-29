b = positionsVector;
x = size(b,2);
plot3(b(1,1:x,1),b(1,1:x,2),b(1,1:x,3),'*');
r = x;
c = zeros(r-1,3);
n = zeros( 1,r-1 );
for i=1:(r-1)
    c(i,1) = b(1,i+1,1) - b(1,i,1);
    c(i,2) = b(1,i+1,2) - b(1,i,2);
    c(i,3) = b(1,i+1,3) - b(1,i,3);
    n(i) = sqrt( c(i,1)^2 + c(i,2)^2 + c(i,3)^2 );
end

nor = [];
figure(2)

for i=1:x
    nor = [ nor sqrt( b(1,i,1)^2 + b(1,i,2)^2 + b(1,i,3)^2 ) ];
end

plot( nor )
