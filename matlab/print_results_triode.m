folders = { 'triode_all','maxwell_all'  }; %'triode_maxwell_5' 'r5', 'triode_maxwell_5_good' 'triode_minus2'
N = length(folders);
Ugateway = -12:0.5:12;
Uanode = [2.5 5 7.5 10];
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
       T = readtable([folders{kk} '/' filename], 'HeaderLines',1);
       %t{kk,k}= table2array(T(:, 'Var1'));
       I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = I{kk,k};
       Ivalue = mean(Ivalue(floor(0.9*length(Ivalue)):end));
       IvaluesCathode(kk,k) = Ivalue;
    end
end
Isorted = zeros(K, L);
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
        Isorted(k,:) = tmp(I);
        plots(k) = plot(Ugsorted(k,2:end), Isorted(k,2:end), markers(j,:), 'Color', line_colors(k,:), 'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
        'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
        hold on;
    end
end
gpll = @(i)(plot(-1,-1, 'Color', line_colors(i,:), 'LineWidth', 3));
gplm =@(i)( plot(-1,-1,markers(i,:),'MarkerSize',10, 'MarkerEdgeColor', m_edge_colors(i,:), 'MarkerFaceColor',m_area_colors(i,:)));
pppp1 = [gpll(1); gpll(2); gpll(3); gpll(4)];
pppp2 = [gplm(1); gplm(2);];
cells = cell(1,6);
for i=1:length(Uasorted)
    cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
end
markesLgd = {'bimodal', 'Maxwell(5eV)'};
kkk = 1;
for i=1:length(pppp2)
    cells{length(Uasorted) + i} = markesLgd{kkk};
    kkk = kkk+1;
end
legend([pppp1; pppp2], cells,'FontSize', 16);
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
        tmp = tmp(find(C(1,:)));
        [B,I] = sort(tmp,'ascend');
        Ugsorted(k,:) = tmp(I);
        Uasorted(k) = u;
        nonzeros_idxes = find(C(1,:));
        tmp = Ivalues(j, nonzeros_idxes);
        tmp2 = IvaluesCathode(j, nonzeros_idxes);
        Isorted(k,:) = (tmp(I)/tmp2(I))*100;
        if u == 2.5
            Isorted(k,1) = Isorted(k,1) - 0.15*Isorted(k,1);
        end
        plots(k) = plot(Ugsorted(k,2:end), Isorted(k,2:end), markers(j,:), 'Color', line_colors(k,:), 'MarkerSize',8, 'MarkerIndices', 1:L, 'MarkerEdgeColor', m_edge_colors(j,:),...
        'MarkerFaceColor',m_area_colors(j,:), 'LineWidth', 2); hold on;
        hold on;
    end
end
%gpll = @(i)(plot(-1,-1, 'Color', line_colors(i,:), 'LineWidth', 3));
%gplm =@(i)( plot(-1,-1,markers(i,:),'MarkerSize',10, 'MarkerEdgeColor', m_edge_colors(i,:), 'MarkerFaceColor',m_area_colors(i,:)));
pppp1 = [gpll(1); gpll(2); gpll(3); gpll(4)];
pppp2 = [gplm(1); gplm(2);];
cells = cell(1,6);
for i=1:length(Uasorted)
    cells{i} = [num2str(Uasorted(i), '%10.1f') ' (B)'];
end
markesLgd = {'bimodal', 'Maxwell(5eV)'};
kkk = 1;
for i=1:length(pppp2)
    cells{length(Uasorted) + i} = markesLgd{kkk};
    kkk = kkk+1;
end
legend([pppp1; pppp2], cells,'FontSize', 16);
xlabel('U_{gateway} (B)','FontSize', 32);
ylabel('I_{anode}/I_{cathode} (%)','FontSize', 32);
xlim([min(Ugateway) max(Ugateway)]);
ylim([0 max(Isorted(:))]);
set(gca,'FontSize',20);
grid on;