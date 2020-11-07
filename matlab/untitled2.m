rootdir = '/home/artoria/results/';
folders = { 'convergence/small_ma', 'convergence/small_bi' };
% 
% M = csvread();
% M2 = csvread();

%figure(1);
%plot_boundaries(tbx, tby );

Files_ma = dir( [rootdir 'convergence/small_ma/'] );
Files_bi = dir( [rootdir 'convergence/small_bi/'] );
expr_cloud = '.*ParticlesCloud.*'; %get collected current file
    
FileNames_ma = {Files_ma.name};
FileNames_bi = {Files_bi.name};
ismatch = @(IN)(~cellfun(@isempty, regexp(IN, expr_cloud, 'match')));
FileNames_clouds_ma = FileNames_ma(ismatch(FileNames_ma));
FileNames_clouds_bi = FileNames_bi(ismatch(FileNames_bi));
K = length(FileNames_clouds_ma);
figures_handled1 = zeros(1,1);
figures_handled2 = zeros(1,1);
for i=1:K
    filename = FileNames_clouds_ma{i};
    filename_to_parse = strsplit(filename, '_');
    gateway_str = filename_to_parse{2};
    anode_str = filename_to_parse{3};
    anode = string(regexp(anode_str,'-?\d?\d.\d\d?','Match'));
    gateway = string(regexp(gateway_str,'-?\d?\d.\d\d?','Match'));
    M1 = csvread([rootdir 'convergence/small_ma/' filename]);
    p1 = sqrt(M1(:,4).^2+M1(:,5).^2+M1(:,6).^2);
    minp1 = min(p1);
    maxp1 = max(p1);
    M2 = csvread([rootdir 'convergence/small_bi/' filename]);
    p2 = sqrt(M2(:,4).^2+M2(:,5).^2+M2(:,6).^2);
    minp2 = min(p2);
    maxp2 = max(p2);
    h1 = zeros(1,2);
    h2 = zeros(1,2);
    q1 = abs(M1(:,7));
    q2 = abs(M2(:,7));
    tiledlayout(2,2);
    h1(1) = nexttile();
    plot_boundaries(tbx, tby ); hold on;
    scatter3( M1(:,1), M1(:,2), p1, 1, p1);
    hcb = colorbar, caxis([ min(minp1,minp2) , max(maxp1,maxp2)]);
    hcb.Label.String = "Magnitude of the generalized impulse";
    hcb.FontSize = 12;
    view(2);
    tt = title( strcat('Maxwell: Generalized impulse: U_{anode}', anode, '(V), U_{gateway}', gateway) );
    h1(2) = nexttile();
    plot_boundaries(tbx, tby ); hold on;
    scatter3( M2(:,1), M2(:,2), p2, 1, p2);
    hcb = colorbar, caxis([ min(minp1,minp2) , max(maxp1,maxp2)]);
    hcb.Label.String = "Magnitude of the generalized impulse";
    hcb.FontSize = 12;
    view(2);
    tt = title( strcat('Bimodal: Generalized impulse: U_{anode}', anode, '(V), U_{gateway}', gateway) );    
    h2(1) = nexttile();
    plot_boundaries(tbx, tby ); hold on;
    scatter3( M1(:,1), M1(:,2), q1, 1, q1);
    hcb = colorbar, caxis([ min( min(q1), min(q2)) , max( max(q1),max(q2))]);
    hcb.Label.String = "Particle's charge abs value (C)";
    hcb.FontSize = 12;
    view(2);
    tt = title( strcat('Maxwell: Absolute value of the charge of particles (C): U_{anode}', anode, '(V), U_{gateway}', gateway) );
    h2(2) = nexttile();
    plot_boundaries(tbx, tby ); hold on;
    scatter3( M2(:,1), M2(:,2), q2, 1, q2);
    hcb = colorbar, caxis([ min( min(q1), min(q2)) , max( max(q1),max(q2))]);
    hcb.Label.String = "Particle's charge abs value (C)";
    hcb.FontSize = 12;
    view(2);
    tt = title( strcat('Bimodal: Absolute value of the charge of particles (C): U_{anode}', anode, '(V), U_{gateway}', gateway) );
    set(h1, 'Colormap', colormap, 'CLim', [ min(minp1,minp2) , max(maxp1,maxp2)])
    set(h2, 'Colormap', colormap, 'CLim', [ min( min(q1), min(q2)) , max( max(q1),max(q2))])
    set(gcf, 'Position', get(0, 'Screensize'));
    exportgraphics(gcf,  strcat('CloudResult/eps/CloudOfParticles_', anode, '_', gateway, '.eps'), 'ContentType', 'image', 'Resolution', 350)
    close(gcf);
end

% x = M(:,1);
% y = M(:,2);
% x2 = M2(:,1);
% y2 = M2(:,2);
% hold on
% p = sqrt(M(:,4).^2+M(:,5).^2+M(:,6).^2);
% p2 = sqrt(M2(:,4).^2+M2(:,5).^2+M2(:,6).^2);
% %gamma = sqrt(ones(size(x))+M(:,4).^2+M(:,5).^2+M(:,6).^2);
% %light_v = 299792458.0;
% %v = light_v .* p ./ gamma;
% q = M(:,7);
% q2 = M2(:,7);
% idx = abs(q) >= 1.60217662e-19;
% idxz = abs(q) == 0;
% cmap = colormap;
% %scatter3( x, y, v, 1, v);
% %colorbar, caxis([ min(v) , max(v)]) % colorbar limits
% scatter3( x, y, p, 3, p);
% colorbar, caxis([ min( min(p), min(p2) ) , max( max(p), max(p2))]) % colorbar limits
% %scatter3( x(idx), y(idx), v(idx), 1, v(idx));
% %colorbar, caxis([ min(v(idx)) , max(v(idx))]) % colorbar limits
% %scatter3(x(idx),y(idx),q(idx), 1, q(idx));
% %plot3(x,y,v, '.');
% %colorbar, caxis([ min(q(idx)) , max(q(idx))]) % colorbar limits
% figure(2);
% plot_boundaries();
% scatter3( x2, y2, p2, 3, p2);
% colorbar, caxis([ min( min(p), min(p2) ) , max( max(p), max(p2))]) % colorbar limits