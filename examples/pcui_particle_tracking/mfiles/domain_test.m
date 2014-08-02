%test depth and slope possibilities
H = 0.15;
h1 = 0.02;
s = 0.2;
Ls = H/s;
Li = (H-h1)/s;

figure;
hold on;
plot([0 0],[-H 0],'k-');
plot([0 Ls],[-H 0],'k-');
plot([0 Ls],[0 0],'k-');
plot([0 Li],[-h1 -h1],'b');
axis image;
axis([-0.1 Ls+0.1 -H-0.1 0.1]);
box on;
hold off;