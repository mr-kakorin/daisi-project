%'modelled_velocity', 'new_nordheim_functions'
folders = { 'newnewnew'};
etalon = [0.05 0.33 1.49 3.25 4.88 11.63 13.1]./400./400;
U = [3.41 4.88 6.44 7.54 9.12 11.19 18.75];
K = length(folders);
L = length(U);
t = cell(K,L);%zeros(1, length(FileNames));
I = cell(K,L);%zeros(1, length(FileNames));
Ua = zeros(K,L);
Ivalues = zeros(K,L);
%'Bimodal', 'Maxwell(0.001 eV)',
lgdForIU = {'calculated'};
lgdForIU{end+1} = 'experimental';
%lgdForDelta = folders;
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
       t{kk,k}= table2array(T(:, 'Var1'));
       I{kk,k} = table2array(T(:, 'Var2'));
       Ivalue = I{kk,k};
       Ivalue = mean(Ivalue(0.9*length(Ivalue):end));
       Ivalues(kk,k) = Ivalue;
       %plot(Ua(kk,k), Ivalue, '.' ,'MarkerSize', 25); hold on;
       %plot(t{kk,k}, I{kk,k}); hold on;
    end
end

[tmp, I] = sort(Ua(1,:),'ascend');
Ua = Ua(:,I);
Ivalues  = Ivalues(:, I);
p = zeros(1, K+1);
for kk=1:K
    %Ivalues(kk,1) = Ivalues(kk,1)+1e-7;
    %Ivalues(kk,2) = Ivalues(kk,2)+6e-7;
    %Ivalues(kk,3) = Ivalues(kk,3)+9e-7;
    %Ivalues(kk,end) = Ivalues(kk,end)+4e-6;
    %Ivalues(kk,end-2) = Ivalues(kk,end-2)-1e-5;
    p(kk) = plot(Ua(kk,:), Ivalues(kk,:), '-o','LineWidth', 2,'MarkerSize', 10); hold on;
end
p(end) = plot(U,etalon, '.' ,'MarkerSize', 35,'LineWidth', 2);
xlabel('U (B)','FontSize', 32);
ylabel('I (A)','FontSize', 32);
lgd = legend(p, lgdForIU,'FontSize', 16);
set(gca,'FontSize',20)
pp = zeros(1,K);
pp2 = zeros(1,K);
figure(2);
for kk=1:K
    %figure(2);
    pp(kk) = plot(Ua(kk,:), (Ivalues(kk,:)-etalon)./etalon * 100 , '-o' ,'MarkerSize', 10,'LineWidth', 2); hold on;
    %figure(3);
    %pp2(kk) = plot(Ua(kk,:), Ivalues(kk,:)-etalon ,'MarkerSize', 25,'LineWidth', 2); hold on;
end
figure(2);
%lgd = legend(pp, lgdForDelta,'FontSize', 16);
xlabel('U (B)','FontSize', 32);
ylabel('\delta (%)','FontSize', 32);
%figure(3);
%xlabel('U (B)','FontSize', 32);
%ylabel('\nabla','FontSize', 32);
%lgd = legend(pp2, lgdForDelta,'FontSize', 16);
set(gca,'FontSize',20)