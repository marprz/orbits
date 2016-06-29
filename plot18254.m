M = out_18254;
v = outvel;
jd0 = 2457023.5;
deltajd = 1/24/12;
mu = 3.98600441800000e+05;
n = size(M,1);
z = [ 0 0 0 ]';
P = [];
for i=1:n
    [ p2, v2, a2 ] = ECEFtoECI( jd0+(i-1)*deltajd, M(i,:)', z, z );
    P(i,:) = p2';
end

x = n/2;
plot3( M(:,1),M(:,2),M(:,3),'.')
hold on;
plot3( M(1:x,1),M(1:x,2),M(1:x,3),'g.')
x = n;
plot3( P(1:x,1),P(1:x,2),P(1:x,3),'r.')

x = n/2;
plot3( P(1:x,1),P(1:x,2),P(1:x,3),'y.')

figure(2)
v_orig = [];
v_calc = [];

n = min(size(M,1),size(v,1));
for i=1:n
    v_orig(i) = norm(v(i,:)/10000);
    r(i) = norm(M(i,:));
    v_calc(i) = sqrt( mu/r(i));
end

plot(v_orig)
hold on;
v_calc
plot( v_calc, 'r')
