folders  = { 'triode_all','maxwell_all' }; %'triode_maxwell_5' 'r5', 'triode_maxwell_5_good' 'triode_minus2'
Ugateway = -12:0.5:12;
Uanode   =  2.5:2.5:10;

rng('default');
markers = {'o-', 's-','^-','v-','>-','<-', '*-', '.-', 'd-', 'x-','p-','h-'};
line_colors   = {'y','m','c','r','g','b','k'};

N = length(folders);
L = length(Ugateway);
K = length(Uanode);

if (N > length(markers))
    return;
end
I  = cell(N,K*L);
Ug = zeros(N, K*L);
Ua = zeros(N, K*L);
U  = zeros(2*N, K*L);

Ivalues        = zeros(N,K*L);
IvaluesCathode = zeros(N,K*L);
IvaluesDiff    = zeros(N,K*L);

indexes = reshape(1:2*N, N, 2)';

for kk=1:N
    Files=dir(folders{kk});
    expressionanode = '.*Collected.*'; %get collected current file
    FileNames = {Files.name};
    ismatch = ~cellfun(@isempty, regexp(FileNames, expressionanode, 'match'));
    FileNamesAnode = FileNames(ismatch);
    for k=1:K*L
       filename  = FileNamesAnode{k}; %extract string from cell
       filename_to_parse = strsplit(filename, '_');
       gateway_str = filename_to_parse{2};
       anode_str = filename_to_parse{3};
       anode = regexp(anode_str,'-?\d?\d.\d\d?','Match');
       gateway = regexp(gateway_str,'-?\d?\d.\d\d?','Match');
       Ua(kk,k) = str2double(anode);
       Ug(kk,k) = str2double(gateway);
       U(indexes(kk,:),k) = [Ua(kk,k) Ug(kk,k)];
       T = readtable([folders{kk} '/' filename], 'HeaderLines',1);
       %I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = table2array(T(:, 'Var2'));
       Ivalue = mean(Ivalue(floor(0.9*length(Ivalue)):end));
       Ivalues(kk,k) = Ivalue;
    end
    expressioncathode = '.*A;.*'; %get flow current file
    ismatch = ~cellfun(@isempty, regexp(FileNames, expressioncathode, 'match'));
    FileNamesCathode = FileNames(ismatch);
    for k=1:K*L
       filename  = FileNamesCathode{k}; %extract string from cell
       filename_to_parse = strsplit(filename, '_');
       gateway_str = filename_to_parse{2};
       anode_str = filename_to_parse{3};
       anode = regexp(anode_str,'\d?\d.\d\d?','Match');
       gateway = regexp(gateway_str,'\d?\d.\d\d?','Match');
       T = readtable([folders{kk} '/' filename], 'HeaderLines',1);
       %I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = table2array(T(:, 'Var2'));
       Ivalue = mean(Ivalue(floor(0.9*length(Ivalue)):end));
       IvaluesCathode(kk,k) = Ivalue;
    end
end
Isorted = zeros(K, L);
Uasorted = zeros(K,1);
Ugsorted = zeros(K,L);
plots = zeros(1, K);

m_edge_colors = repelem([0,0,0], N, 1, 1);
m_area_colors = rand(N, 3);

f1 = figure(1);
for j=1:N
    for k=1:K
        u = Uanode(k);
        C = U(indexes(j,:), :) == u;
        tmp = U(2,:);
        tmp = tmp(C(1,:));
        [B,I] = sort(tmp,'ascend');
        Ugsorted(k,:) = tmp(I);
        Uasorted(k) = u;
        nonzeros_idxes = find(C(1,:));
        tmp = Ivalues(j, nonzeros_idxes);
        Isorted(k,:) = tmp(I);
        plots(k) = plot(Ugsorted(k,:), Isorted(k,:), markers{j}, 'Color', line_colors{k},...
            'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
            'MarkerFaceColor', m_area_colors(j,:), 'LineWidth', 2);
        hold on;
    end
end
gpll = @(i)(plot(-1,-1, 'Color', line_colors{i}, 'LineWidth', 3));
gplm = @(i)( plot(-1,-1, markers{i}, 'MarkerSize', 10, 'MarkerEdgeColor', m_edge_colors(i,:), 'MarkerFaceColor',m_area_colors(i,:)));
pppp1 = arrayfun(gpll, 1:K);
pppp2 = arrayfun(gplm, 1:N);
cells = cell(1, K+N);
for i=1:length(Uasorted)
    cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
end
kkk = 1;
for i=1:length(pppp2)
    cells{length(Uasorted) + i} = folders{kkk};
    kkk = kkk+1;
end
legend([pppp1 pppp2], cells, 'FontSize', 16);

xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{anode} (A)','FontSize', 32);
xlim([min(Ugateway) max(Ugateway)]);
ylim([0 max(Ivalues(:))]);
set(gca,'FontSize',20);
grid on;

figure(2);
for j=1:N
    for k=1:K
        u = Uanode(k);
        C = U(indexes(j,:), :) == u;
        tmp = U(2,:);
        tmp = tmp(C(1,:));
        [B,I] = sort(tmp,'ascend');
        Ugsorted(k,:) = tmp(I);
        Uasorted(k) = u;
        nonzeros_idxes = find(C(1,:));
        tmp = Ivalues(j, nonzeros_idxes);
        tmp2 = IvaluesCathode(j, nonzeros_idxes);
        Isorted(k,:) = (tmp(I)/tmp2(I))*100;
        plots(k) = plot(Ugsorted(k,:), Isorted(k,:), markers{j}, 'Color', line_colors{k},...
            'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
            'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
        hold on;
    end
end
pppp1 = arrayfun(gpll, 1:K);
pppp2 = arrayfun(gplm, 1:N);
legend([pppp1 pppp2], cells,'FontSize', 16);
xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{anode}/I_{cathode} (%)','FontSize', 32);
xlim([min(Ugateway) max(Ugateway)]);
ylim([0 max(Isorted(:))]);
set(gca,'FontSize',20);
grid on;