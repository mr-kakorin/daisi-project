rootdir = '/home/artoria/results/';
%folders = { 'convergence/small_ma', 'convergence/big_ma', 'convergence/small_bi', 'convergence/big_bi' };
folders = { 'convergence/small_ma', 'convergence/small_bi' };%'triode_maxwell_5' 'r5', 'triode_maxwell_5_good' 'triode_minus2'
%markesLgd = {'small-maxwell', 'big-maxwell', 'small-bimodal', 'big-bimodal'};      
markesLgd = {'maxwell 8.2792 Ev', 'bimodal exprmt'};    
N = length(folders);
Ugateway = 1:20;
Uanode = [9];
L = length(Ugateway);
K = length(Uanode);
t = cell(N,K*L);
I = cell(N,K*L);
Ug = zeros(N,K*L);
Ua = zeros(N,K*L);
U = zeros(2*N, K*L);
Ivalues = zeros(N,K*L);
IvaluesCathode = zeros(N,K*L);
IvaluesDiff = zeros(N,K*L);
indexes=[1 2; 3 4; 5 6; 7 8; 9 10];

for kk=1:N
    folders_kk = [ rootdir folders{kk} ];
    Files=dir( folders_kk );
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
       T = readtable([ folders_kk '/' filename], 'HeaderLines',1);
       %t{kk,k}= table2array(T(:, 'Var1'));
       I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = I{kk,k};
       Ivalue = mean(Ivalue(floor(0.9*length(Ivalue)):end));
       Ivalues(kk,k) = Ivalue;
    end
    expressioncathode = '.*A;.*'; %get collected current file
    ismatch = ~cellfun(@isempty, regexp(FileNames, expressioncathode, 'match'));
    FileNamesCathode = FileNames(ismatch);
    for k=1:K*L
       filename  = FileNamesCathode{k}; %extract string from cell
       filename_to_parse = strsplit(filename, '_');
       gateway_str = filename_to_parse{2};
       anode_str = filename_to_parse{3};
       anode = regexp(anode_str,'\d?\d.\d\d?','Match');
       gateway = regexp(gateway_str,'\d?\d.\d\d?','Match');
       %Ua(kk,k) = str2double(anode);
       %Ug(kk,k) = str2double(gateway);
       %U(indexes(kk,:),k) = [Ua(kk,k) Ug(kk,k)];
       T = readtable([folders_kk '/' filename], 'HeaderLines',1);
       %t{kk,k}= table2array(T(:, 'Var1'));
       I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = I{kk,k};
       Ivalue = mean(Ivalue(floor(0.9*length(Ivalue)):end));
       IvaluesCathode(kk,k) = Ivalue;
    end
end
Isorted = zeros(N, K, L);
Uasorted = zeros(K,1);
Ugsorted = zeros(K,L);
plots = zeros(1, K);
markers = ['-s'; '-o'; '-^'; '-v'; '-d'];
m_edge_colors = [[0,0,0];...
                 [0,0,0];...
                 [0,0,0];...
                 [0,0,0];...
                 [0,0,0]];
m_area_colors = [[1,0,0];...
                 [0,1,0];...
                 [0,0,1];...
                 [0.5, 0, 0.5];...
                 [0, 0.25, 0.5]];
line_colors = [[0, 0.4470, 0.7410];...
               [0.4940, 0.1840, 0.5560];...
               [0.6350, 0.0780, 0.1840];...
               [0.9290, 0.6940, 0.1250]];
         
f1 = figure(1);
for j=1:N
    for k=1:K
        u = Uanode(k);
        C = U(indexes(j,:), :) == u;
        tmp = U(2,:);
        tmp = tmp(find(C(1,:)));
        [B,I] = sort(tmp,'ascend');
        Ugsorted(k,:) = tmp(I);
        Uasorted(k) = u;
        nonzeros_idxes = find(C(1,:));
        tmp = Ivalues(j, nonzeros_idxes);
        Isorted(j, k,:) = tmp(I);
        plots(k) = plot(Ugsorted(k,1:end), squeeze( Isorted(j, k,1:end) ), markers(j,:), 'Color', line_colors(k,:), 'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
        'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
        hold on;
    end
end
pppp1 = zeros( size( Uanode ) );
pppp2 = zeros( size( folders ) );

gpll = @(i)(plot(-1,-1, 'Color', line_colors(i,:), 'LineWidth', 3));
gplm = @(i)( plot(-1,-1,markers(i,:),'MarkerSize',10, 'MarkerEdgeColor', m_edge_colors(i,:), 'MarkerFaceColor',m_area_colors(i,:)));
for i=1:length(pppp1)
   pppp1(i) = gpll(i);
end
for i=1:length(pppp2)
   pppp2(i) = gplm(i);
end

cells = cell(1,N+K);
for i=1:length(Uasorted)
    cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
end
kkk = 1;
for i=1:length(pppp2)
    cells{length(Uasorted) + i} = markesLgd{kkk};
    kkk = kkk+1;
end
legend([pppp1 pppp2], cells,'FontSize', 16);
xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{anode} (A)','FontSize', 32);
xlim([0 max(Ugateway)]);
ylim([0 max(Ivalues(:))]);
set(gca,'FontSize',20);
grid on;


figure(2);
IsortedCathode = zeros(N, K, L);
for j=1:N
    for k=1:K
        u = Uanode(k);
        C = U(indexes(j,:), :) == u;
        tmp = U(2,:);
        tmp = tmp(find(C(1,:)));
        [B,I] = sort(tmp,'ascend');
        Ugsorted(k,:) = tmp(I);
        Uasorted(k) = u;
        nonzeros_idxes = find(C(1,:));
        tmp = IvaluesCathode(j, nonzeros_idxes);
        IsortedCathode(j, k, :) = tmp(I);
        plots(k) = plot(Ugsorted(k,1:end), squeeze(IsortedCathode(j, k, 1:end)), markers(j,:), 'Color', line_colors(k,:), 'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
        'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
        hold on;
    end
end
for i=1:length(pppp1)
   pppp1(i) = gpll(i);
end
for i=1:length(pppp2)
   pppp2(i) = gplm(i);
end

for i=1:length(Uasorted)
    cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
end
kkk = 1;
for i=1:length(pppp2)
    cells{length(Uasorted) + i} = markesLgd{kkk};
    kkk = kkk+1;
end
legend([pppp1 pppp2], cells,'FontSize', 16);
xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{cathode} (A)','FontSize', 32);
xlim([0 max(Ugateway)]);
ylim([0 max(IvaluesCathode(:))]);
set(gca,'FontSize',20);
grid on;

figure(3);

pppp = zeros(1,2);
pppp(1) = plot( Ugsorted(1, :), abs( squeeze( Isorted(1, 1, :) )./squeeze(Isorted(2,1,:))*100), 'LineWidth', 2 );
hold on;
pppp(2) = plot( Ugsorted(1, :), abs( squeeze(IsortedCathode(1,1, :))./squeeze(IsortedCathode(2,1,:))*100), 'LineWidth', 2 );
grid on;

borderx = [ 0 max( max(abs( squeeze( Isorted(1, 1, :) )./squeeze(Isorted(2,1,:))*100)), max(abs( squeeze(IsortedCathode(1,1, :))./squeeze(IsortedCathode(2,1,:))*100)) )];

xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{maxwell}/I_{bimodal} (%)','FontSize', 32);
xlim([0 max(Ugateway)]);
ylim(borderx);
legend( pppp, {'anode', 'cathode'},'FontSize', 16 );

% figure(4);
% for j=1:N
%     for k=1:K
%         u = Uanode(k);
%         C = U(indexes(j,:), :) == u;
%         tmp = U(2,:);
%         tmp = tmp(find(C(1,:)));
%         [B,I] = sort(tmp,'ascend');
%         Ugsorted(k,:) = tmp(I);
%         Uasorted(k) = u;
%         nonzeros_idxes = find(C(1,:));
%         tmp = Ivalues(j, nonzeros_idxes);
%         tmp2 = IvaluesCathode(j, nonzeros_idxes);
%         Isorted(k,:) = (tmp(I)/tmp2(I))*100;
%         plots(k) = plot(Ugsorted(k,2:end), Isorted(k,2:end), markers(j,:), 'Color', line_colors(k,:), 'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
%         'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
%         hold on;
%     end
% end
% %gpll = @(i)(plot(-1,-1, 'Color', line_colors(i,:), 'LineWidth', 3));
% %gplm =@(i)( plot(-1,-1,markers(i,:),'MarkerSize',10, 'MarkerEdgeColor', m_edge_colors(i,:), 'MarkerFaceColor',m_area_colors(i,:)));
% for i=1:length(pppp1)
%    pppp1(i) = gpll(i);
% end
% for i=1:length(pppp2)
%    pppp2(i) = gplm(i);
% end
% 
% for i=1:length(Uasorted)
%     cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
% end
% %markesLgd = {'bimodal', 'Maxwell 8.292eV', 'bimodal2', 'Maxwell2 8.292eV'};
% kkk = 1;
% for i=1:length(pppp2)
%     cells{length(Uasorted) + i} = markesLgd{kkk};   
%     kkk = kkk+1;
% end
% legend([pppp1 pppp2], cells,'FontSize', 16);
% xlabel('U_{gateway} (B)','FontSize', 32);
% ylabel('I_{anode}/I_{cathode} (%)','FontSize', 32);
% xlim([min(Ugateway) max(Ugateway)]);
% ylim([0 max(abs(Isorted(:)))+max(abs(Isorted(:)))*0.1]);
% set(gca,'FontSize',20);
% grid on;