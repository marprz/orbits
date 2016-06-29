a = out_file;

for i=1:size(a,1)
    plot3(a(i,:,1),a(i,:,2),a(i,:,3));
    hold on
    plot3(a(i,1,1),a(i,1,2),a(i,1,3),'*r');
    plot3(a(i,288,1),a(i,288,2),a(i,288,3),'*k');
end
