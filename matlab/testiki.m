tiledlayout(2,2);
h1 = zeros(1,2);
h2 = zeros(1,2);
q1 = abs(M1(:,7));
q2 = abs(M2(:,7));
h1(1) = nexttile();
plot_boundaries(tbx, tby ); hold on;
scatter3( M1(:,1), M1(:,2), p1, 1, p1);
hcb = colorbar, caxis([ min(minp1,minp2) , max(maxp1,maxp2)]);
hcb.Label.String = "Magnitude of the generalized impulse";
title('Generalized impulse')
view(2);
h1(2) = nexttile();
plot_boundaries(tbx, tby ); hold on;
scatter3( M2(:,1), M2(:,2), p2, 1, p2);
hcb = colorbar, caxis([ min(minp1,minp2) , max(maxp1,maxp2)]);
hcb.Label.String = "Magnitude of the generalized impulse";
view(2);
h2(1) = nexttile();
plot_boundaries(tbx, tby ); hold on;
scatter3( M1(:,1), M1(:,2), q1, 1, q1);
hcb = colorbar, caxis([ min( min(q1), min(q2)) , max( max(q1),max(q2))]);
hcb.Label.String = "Particle charge abs value (Кл)";
view(2);
h2(2) = nexttile();
plot_boundaries(tbx, tby ); hold on;
scatter3( M2(:,1), M2(:,2), q2, 1, q2);
hcb = colorbar, caxis([ min( min(q1), min(q2)) , max( max(q1),max(q2))]);
hcb.Label.String = "Particle charge abs value (Кл)";
view(2);

set(h1, 'Colormap', colormap, 'CLim', [ min(minp1,minp2) , max(maxp1,maxp2)])
set(h2, 'Colormap', colormap, 'CLim', [ min( min(q1), min(q2)) , max( max(q1),max(q2))])
set(gcf, 'Position', get(0, 'Screensize'));