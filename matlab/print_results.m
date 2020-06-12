folders = { 'd7', 'd8' };

etalon = [0.05 0.33 1.49 3.25 4.88 11.63 13.1]./400./400;
U = [3.41 4.88 6.44 7.54 9.12 11.19 18.75];

rng('default');
markers = {'o-', 's-','^-','v-','>-','<-', '*-', '.-', 'd-', 'x-','p-','h-'};
line_colors   = {'y','m','c','r','g','b','k'};
K = length(folders);
L = length(U);

if (K > length(markers))
    return;
end

I = cell(K,L);
Ua = zeros(K,L);
Ivalues = zeros(K,L);
for kk=1:K
    %figure;
    Files=dir(folders{kk});
    expression = '.*A;.*';
    FileNames = {Files.name};
    ismatch = ~cellfun(@isempty, regexp(FileNames, expression, 'match'));
    FileNames = FileNames(ismatch);
    for k=1:L
       filename  = FileNames{k}; %extract string from cell
       filename_to_parse = strsplit(filename, '_');
       filename_to_parse = filename_to_parse{2};
       anode = regexp(filename_to_parse,'\d?\d.\d\d?','Match');
       Ua(kk,k) = str2double(anode);
       T = readtable([folders{kk} '/' filename], 'HeaderLines',1);
       I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = I{kk,k};
       Ivalue = mean(Ivalue(0.9*length(Ivalue):end));
       Ivalues(kk,k) = Ivalue;
    end
end

[tmp, I] = sort(Ua(1,:),'ascend');
Ua = Ua(:,I);
Ivalues  = Ivalues(:, I);
p = zeros(1, K+1);
for kk=1:K
    p(kk) = plot(Ua(kk,:), Ivalues(kk,:), markers{kk},'LineWidth', 2,'MarkerSize', 10); hold on;
end
p(end) = plot(U, etalon, markers{kk+1} ,'MarkerSize', 10, 'LineWidth', 2);
xlabel('U (B)','FontSize', 32);
ylabel('I (A)','FontSize', 32);
legend(p, [folders "etalon"],'FontSize', 16);
set(gca,'FontSize',20)

figure(2);
pp = zeros(1,K);
for kk=1:K
    pp(kk) = plot(Ua(kk,:), (Ivalues(kk,:)-etalon)./etalon * 100 , markers{kk} ,'MarkerSize', 10,'LineWidth', 2); hold on;
end
legend(pp, folders,'FontSize', 16);
xlabel('U (B)','FontSize', 32);
ylabel('\delta (%)','FontSize', 32);
set(gca,'FontSize',20)