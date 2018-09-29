surf(X(1:9,3:19),Y(1:9,3:19),Zl(1:9,3:19)), colorbar;
xlabel('log(ratio of db1/db2)');
ylabel('thickness');
zlabel('Percent Error in db2');
set(gca,'YDir','reverse');
c = colorbar;
caxis([-20 20]);  