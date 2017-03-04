function [S2min]=model9p(r_in,d_in,p_in,su1_in,su2_in)
% F4 |(1,3,1),(2,3,4),(1,3,1),(2,3,4),(1,3,1),(2,3,4,5) : (1-2-(1-2)-1-2)

flag1 = 1;
if flag1==0
    clc
    clear
    
    load('db_data2.mat')
    load('id_data2_1.mat')
    load('id_data2_2.mat')
    load('p_data2.mat');
    load('su1_data2.mat');
    load('su2_data2.mat');
    load('r_data2.mat');
    load('d_data2.mat')
end

%n=312; %10
[~,n]=size(p_in);
n = min(n,300);
n_glob = 1;

n=ceil(n/n_glob);

flag_branch = 0;
flag_mneh = 0; %mneh
flag_mneh_fl = 0; %mneh-fl
flag_1212 = 1; %1212-fl
flag_1212_fl = 0; %1212-fl
Kross_fl = 0.001;
flag_SACO1 = 0; %SACO полная оптимизация
flag_SACO_fl_1 = 0; %SACO-fl
flag_KUMZ = 0; %SACO-fl
flag_SACO2 = 0; %SACO частичная оптимизация по 12-12
flag_fl3_2 = 0; %SACO-fl
flag_GA1 = 0; %ga полная оптимизация
flag_GA2 = 0; %ga частичная оптимизация по 12-12

flag_full = 1; %тележки и выход

flag_131 = 0;
flag_jhons =0;
H=10000;

[~,n_max]=size(p_in); %312
m=5;

%global record;
%global xmin;
%global ccount; %счетчик заходов для теста

% global F_max;
% global a_max;

%record = inf;
%xmin=[];
%ccount=0;
Smin = zeros(n_glob,n*19);

for glob_i=1:n_glob
    
    %p = ceil(rand(19,5)*10);        %пример 1
    %p=p_in;
    % p = [54  53    52  50  50  76  50  51  51  53;  %1 - 1
    %     9     9      9     0    0     0   0     0    0    9;      %3  2
    %     27   26   26   0    0     0   0     0    0    53;   %1   3
    %     17   17   34   36 36   34  35  35  35  35;  %2   4
    %     9     9     9      9   9     9    9     9    9    9;   %3     5
    %     6     5     5      5   4     0    0     0    0    5;     %4   6
    %     28   26   26  50  49  76  49   50  50  53;   %1   7
    %     9     9      9    0    0     0     0    0    0    9;     %3   8
    %     58   54   54  0    0     0     0    0    0    55;     %1 9
    %     17   17   34  36  36  35  35   35  35  35;  %2   10
    %     9     9     9     9    9     9    9     9    9    9;   %3     11
    %     0     0     0     0    0     0    0     0    5    0;   %4     12
    %     63   56  56   52  50   51  33  33  52  0;   %1     13
    %     9     9     9     9    9      9   9     9    9    0;      %3  14
    %     0     0     0     18  17   52  33  34  0    0;     %1   15
    %     17   17  34    36 36   34  38  35  35  0;  %2      16
    %     9     9     9      9   9     9     9    9    9    0;   %3     17
    %     33   0    30    0    0     0     0    0    0    0;   %4    18
    %     0     41   0     48   39  36  35  38  46  40];   %5 19
    
    %снизим размерность иначе долго считает
    %p=round(p_in(:,1:n)*100)/100;
    
    p=p_in(:,1+(glob_i-1)*n:n*glob_i); %1:n
    
    su=cell(1,19);
    
    % проверяем открыт ли пул работ
    %if isempty(gcp('nocreate')) %нет
    %    parpool('local',4); % загрузим 4 процессора
    %end
    
    %инициализация
    %ast1= tic;
    for ii=1:19 %parfor % бессмысленнен
        su{ii}=zeros(n,n); %zeros(n_max,n_max)
    end
    %aend1 = toc(ast1);
    
    % валки
    su{1}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{3}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{7}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{9}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{13}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{15}=su1_in(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    % ножи
    su{6}=su2_in{1}(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{12}=su2_in{2}(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    su{18}=su2_in{3}(1+(glob_i-1)*n:n*glob_i,1+(glob_i-1)*n:n*glob_i);
    
    %r=r_in(1+(glob_i-1)*n:n*glob_i);
    %пока не уверен в равильности выгрузки r
    r = zeros(1,n); %release dates - добавлены в алгоритм
    
    %r = ceil(rand(1,n)*2300);
    %r = [150 0 195 214 156 174 171 90 150 39];
    
    %тест
    %r(2) = 600;
    %r(5) = 600;
    
    %d = ones(1,n)*4700; %due dates
    %d = ceil(rand(1,n)*4700); %4700 9300
    %d = [3082 168 3991 4390 3191 3562 3493 1844 3081 805];
    %d = [5431 14387 14321 7229 10789 1442 1026 10086 14805 17747];
    %d = [3881 462 8396 8787 4566 4551 3141 8371 3434 1035 7257 3625 2248 3757 898 1228 8762 8893 5350 556 2184 3285 7638 144 401 1572 6037 6806 6025 4194];
    
    d = d_in(1+(glob_i-1)*n:n*glob_i);
    %d= max(d,0);
    
    ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5]; % последовательность РЦ
    mq =[1 473 1 18 473 1 1 473 1 18 473 1 1 473 1 18 473 1 2]; % кол-во машин в РЦ
    mq_u = [1 18 473 1 2];
    
    %начнем с примитивных алгоритмов
    % 1. EDD
    [~, I] = sort(d);
    %[~, Ir] = sort(r);
    
    % важно! KumzFullS распараллелен
    S_EDD = KumzFullS(I,ma); % (perm) расписание - послед-й проход заказов:
    %S_EDD = KumzFullSp(I,ma); % (perm) расписание - послед-й проход заказов:
    % важно! XieTargFun7/12f не параллелится
    %[T_EDD1_1,~,~, end_time_ex_EDD1_1]= XieTargFun7(p, S, n, m, ma, d, su);
    %добавлена множественность РЦ
    %[T_EDD1,~,~, end_time_ex_EDD1]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
    %добавлена проверка r
    %tmp
    %[T_EDD1,~,~, end_time_ex_EDD1]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
    [T_EDD,~,~, end_time_ex_EDD]= KUMZTargFun19f(p, S_EDD, n, m, ma, d, su, mq, r, mq_u);
    Cmax_EDD = max(max(end_time_ex_EDD));
    %[start_time_EDD1_full, end_time_EDD1_full]= KUMZ_data_load3(p,S, n, m, ma, d, su, mq, r, mq_u, db_data, id_in);
    %полное расписание строится только для оптимальной функции!
    %[end_time_ex_EDD_full]= KUMZ_data_load(p,S, n, m, ma, d, su, mq, db_data,id_in);
    
    % важно! выгоднее проходить последовательно по каждому агрегату
    %S = KumzFullS2(I,ma); % (perm) расписание - послед-й на каждом рц:
    %[T_EDD2_1,~,~, end_time_ex_EDD2_1] =  XieTargFun7(p, S, n, m, ma, d, su);
    %[T_EDD2_2,~,~, end_time_ex_EDD2_2] = XieTargFun12f(p, S, n, m, ma, d, su, mq);
    %[T_EDD2_2,~,~, end_time_ex_EDD2_2] = KUMZTargFun16f(p, S, n, m, ma, d, su, mq, r);
    
    % 2. максимальная длительность
    p1=sum(p(:,1:n));
    [~,I]=sort(p1,'ascend');
    S_maxp = KumzFullS(I,ma);
    %S_maxp = KumzFullSp(I,ma);
    %S = KumzFullS2(I,ma); % (perm) расписание - послед-й на каждом рц:
    %[T_maxp,~,~, end_time_ex_maxp] = XieTargFun7(p, S, n, m, ma, d, su);
    %[T_maxp,~,~, end_time_ex_maxp]=XieTargFun12f(p, S, n, m, ma, d, su, mq);
    %[T_maxp,~,~, end_time_ex_maxp] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
    [T_maxp,~,~, end_time_ex_maxp] = KUMZTargFun19f(p, S_maxp, n, m, ma, d, su, mq, r, mq_u);
    k1 = 0.3*max(max(end_time_ex_maxp))/T_maxp; %необходимо шкалирование
    FT_Maxp = (T_maxp - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_maxp)) - Cmax_EDD);
    
    %3. пробуем решить SACO полностью - perm расписание
    if flag_SACO1==1
        S121_1 = [];
        S121_2 = [];
        [~,S_opt_SACO,~,~,~]=SACO_KUMZ_12_x_12_4p(p, n, m, ma, d, su,S121_1,S121_2, r,flag_SACO_fl_1,flag_KUMZ, mq, mq_u, T_EDD, Cmax_EDD,k1);
        S_SACO = KumzFullS(S_opt_SACO,ma);
        
        [T_SACO,~,~, end_time_ex_SACO] = KUMZTargFun19f(p, S_SACO, n, m, ma, d, su, mq, r, mq_u);
        FT_SACO = (T_SACO - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_SACO)) - Cmax_EDD);
    else
        FT_SACO = inf;
        S_SACO = S_EDD;
    end
    
    %3. пробуем решить GA полностью - perm расписание
    if flag_GA1 == 1
        %важно! стандартный ГА не подходит, т.к. не имеет параметров отвечающих
        %за различность - для примера использования оставим
        %         FitnessFunction = @(S) KUMZTargFun19fga(p, S, n, m, ma, d, su, mq, r, mq_u);
        %         NonLinConstr = @(S) GA_unique_set2(S,n);
        %
        %         [x,fval] = ga(FitnessFunction,n*19,[],[],[],[],zeros(1,n*19),ones(1,n*19),[],1:n*19); %1:n %ones(n),ones(1,n)*sum(1:n) %1:n*19
        
        S121_1 = [];
        S121_2 = [];
        [ ~, S_opt_GA ] = HGA_KUMZ( p, n, m, ma, d, su,S121_1,S121_2, r,mq, mq_u, flag_KUMZ,T_EDD, Cmax_EDD,k1);
        S_GA = KumzFullS(S_opt_GA,ma);
        
        [T_GA,~,~, end_time_ex_GA] = KUMZTargFun19f(p, S_GA, n, m, ma, d, su, mq, r, mq_u);
        FT_GA = (T_GA - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_GA)) - Cmax_EDD);
    else
        FT_GA = inf;
        S_GA = S_EDD;
    end
    
    % 5. MNEH - FL целиком для задачи
    if flag_mneh ==1
        d_m = d;
        [~, I] = sort(d_m);
        
        S_MEDD = KumzFullS2(I,ma);
        
        S0 = S_MEDD;
        nk = numel(S_MEDD);
        %S=[];
        S = zeros(1,nk);
        %S_new_par = [];
        %ri=1;
        %while ri<=2*n
        %ast1=tic;
        for  ri = 1:nk   %этот цикл не параллелится!
            sr = S0(ri);
            sr2 = S0(ri+1:numel(S0));
            %T_min = inf;
            %S_min = [];
            %C_min = inf;
            T_S_new=ones(1,ri)*inf; %numel(S)
            SS1 = S(1:ri-1);
            parfor l=1:ri %parfor эффективен  % l=1:numel(S) +1  % +1
                if l==1
                    S_new_par = [sr SS1];
                elseif l == ri  %numel(S) % +1
                    S_new_par = [SS1 sr];
                else
                    S_new_par = [SS1(1:l-1) sr SS1(l:ri-1)]; %S(l:numel(S))
                end
                %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
                S_new = [S_new_par sr2];
                %S_new_full = JeongFullS(S_new);
                %S_new_full = JeongFullSp(S_new);
                %важно! вместо d - d1
                %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
                % важно! добавил r
                [T_S_new(l),~,~,~]= KUMZTargFun19f(p, S_new, n, m, ma, d, su, mq, r, mq_u); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
            end
            [T_min, ind_m] = min(T_S_new);
            T_min_old = T_min;
            %S=S_min;
            if ind_m == 1
                %S = [sr S];
                S(2:ri) = S(1:ri-1);
                S(1) = sr;
            elseif ind_m == ri %numel(S) + 1
                %S = [S sr];
                S(ri) = sr;
            else
                %S = [S(1:ind_m-1) sr S(ind_m:numel(S))];
                S(ind_m+1:ri)=S(ind_m:ri-1);
                S(ind_m) = sr;
            end
            
            if flag_mneh_fl ==1
                %конец MNEH начало FL
                if ri>=3
                    %S_min = [];
                    T_S_new=ones(ri-1,ri)*inf; %numel(S)
                    SS1 = S(1:ri);
                    parfor ii=1:ri-1 %эффективен parfor - полный перебор всех соседей
                    %for ii=1:1 %перебор только 1-го
                        T_S_new(ii,:) = S_xchange19(ii,SS1,ri,sr2, p, n, m, ma, d, su, r, mq, mq_u);
                    end
                    [T_S_new1,ind1] = min(T_S_new);
                    [T_min,ind2] = min(T_S_new1);
                    i=ind1(ind2);
                    j=ind2;
                    if T_min<T_min_old
                        tempo = S(j);
                        S(j)=S(i);
                        S(i) = tempo;
                    end
                else
                end
            end
        end
        
        S_MNEH = S;
        [T_MNEH,~,~, end_time_ex_MNEH] = KUMZTargFun19f(p, S_MNEH, n, m, ma, d, su, mq, r, mq_u);
        FT_MNEH = (T_MNEH - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_MNEH)) - Cmax_EDD);
    else
        FT_MNEH = inf;
        S_MNEH = S_EDD;
    end
    
    %--------------------------------------------------------------------------
    % 6. алгоритм NEH
    S=1:n;
    % важно! NawazHehCmax не параллелится
    %[~,NEH3]=NawazHehCmax(p',m,ma,S);
    %добавлена проверка r
    %важно! в итоге результат еще хуже за счет того, что не исп-ся внутрянняя
    %оптимизация TargFun15f, пока оставим так
    [~,NEH3]=NawazHehCmax2(p',m,ma,S,r);
    S_NEH_f = KumzFullS(NEH3,ma); % (perm) расписание - послед-й проход заказов:
    %S_NEH_f = KumzFullSp(NEH3,ma); % parfor неэффективен
    %S = KumzFullS2(NEH3,ma); % (perm) расписание - послед-й на каждом рц:
    %[T_NEH_f,~,~, end_time_ex_NEH_f]= XieTargFun7(p, S, n, m, ma, d, su); %в1
    %[T_NEH_f,~,~, end_time_ex_NEH_f]= XieTargFun12f(p, S, n, m, ma, d, su,mq); %в2
    %[T_NEH_f,~,~, end_time_ex_NEH_f] = KUMZTargFun17f(p, S, n, m, ma, d, su,mq, r, mq_u); %%в3
    [T_NEH_f,~,~, end_time_ex_NEH_f] = KUMZTargFun19f(p, S_NEH_f, n, m, ma, d, su, mq, r, mq_u);
    FT_NEH_f = (T_NEH_f - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_NEH_f)) - Cmax_EDD);
    
    %--------------------------------------------------------------------------
    % пробуем решить как хотели: (1-2-1-2)-1
    if flag_1212 == 1
        % в. 1. просто сложив (1-3-1)=1,(2,3,4)=2,(1-3-1)-2=1,(2,3,4)-2=2,
        % ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5];
        %---------------------------------------------------------
        %
        if flag_131 == 1
            % пробуем решить задачу (1,3,1)-1
            a1=p(1,:);
            b1=p(2,:);
            c1=p(3,:);
            
            a=sum(a1);
            b=sum(b1);
            c=sum(c1);
            
            pj=[a1;b1]';
            S1= jhons(pj);
            Cab = WangTargFun(a1,b1,c1,S1);
            if Cab == a + c || flag_jhons ==1
                xmin1 = S1;
                %record1 = Cab;
                %exit=0;
            else
                [F,xmin1]=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
                %xmin1 = xmin;
                %record1 = record;
            end
            
            % пробуем решить задачу (1,3,1)-2
            a1=p(7,:);
            b1=p(8,:);
            c1=p(9,:);
            
            a=sum(a1);
            b=sum(b1);
            c=sum(c1);
            
            pj=[a1;b1]';
            S1= jhons(pj);
            Cab = WangTargFun(a1,b1,c1,S1);
            if Cab == a + c || flag_jhons ==1
                xmin2 = S1;
                %record2 = Cab;
                %exit=0;
            else
                [F,xmin2]=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
                %xmin2 = xmin;
                %record2 = record;
            end
            % используя решения задач (1,3,1)-1 и (1,3,1)-2
            % модифицируем d
            [~,I]=sort(xmin1);
            d1=d(xmin1);
            
            
            for ii=1:n
                d1(ii)=d1(ii)-H/ii+H/n;
            end
            d2=d1(I);
            
            % модифицируем d по второй тройке
            [~,I]=sort(xmin2);
            d1=d2(xmin2);
            %H=200;
            
            for ii=1:n
                d1(ii)=d1(ii)-H/ii+H/n;
            end
            d2=d1(I);
        else
            d2=d;
        end
        
        p1=zeros(4,n);
        p1(1,:) = sum(p(1:3,:));
        p1(2,:) = sum(p(4:6,:));
        p1(3,:) = sum(p(7:9,:));
        p1(4,:) = sum(p(10:12,:));
        %
        su1 = cell(1,4);
        %for ii=1:4 %parfor
        su1{1}(:,:)=su{1} + su{2} + su{3}; %sum(su{ii}(1:3,:));
        su1{2}(:,:)=su{4} + su{5} + su{6}; %sum(su{ii}(4:6,:));
        su1{3}(:,:)=su{7} + su{8} + su{9}; %sum(su{ii}(7:9,:));
        su1{4}(:,:)=su{10} + su{11} + su{12}; %sum(su{ii}(10:12,:));
        %end
        
        %d1 = d - sum(p(13:19,:));
        %т.к. d уже с учетом реальных длительностей
        d1 = d2 - p(13,:) - p(14,:)*473 - p(15,:) - p(16,:)*18 - p(17,:)*473 - p(18,:) - p(14,:)*2; %d
        
        %применим алгоритм NEH для начального приближения
        %важно! так себе приближение
        S_NEH=1:n;
        m1=2;
        ma1=[1 2 1 2];
        %для NEH необх транспонировать
        %[~,NEH1]=NawazHehCmax(p1',m1,ma1,S_NEH);
        [~,NEH]=NawazHehCmax2(p1',m1,ma1,S_NEH,r);
        
        %S_NEH = [];
        S_NEH = zeros(1,2*n);
        %ast1=tic;
        for ii=1:n %parfor бесмысленнен
            S_NEH(2*ii-1) = NEH(ii);
            S_NEH(2*ii) = NEH(ii);
            %tmp = [NEH(ii) NEH(ii)];
            %S_NEH = [S_NEH tmp]; %ускорение
        end
        %aend1=toc(ast1);
        
        %ast1=tic;
        S = JeongFullS(S_NEH);
        %S = JeongFullSp(S_NEH); %паралл версия бесмысленна
        %aend1=toc(ast1);
        %[T_NEH_1212,~,~,~]= XieTargFun7(p1, S, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные
        % важно! добавил r
        [T_NEH_1212,~,~,~]= XieTargFun7r(p1, S, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные
        
        xmin=S_NEH;
        record = T_NEH_1212; %еще как используется
        
        %применим алгоритм модифицированный FL - MFL для начального приближения
        %MEDD
        d_m = d1 - (p1(3,:) + p1(4,:));
        [~, I] = sort(d_m);
        
        %S_MEDD = [];
        S_MEDD = zeros(1,2*n);
        %ast1=tic;
        for ii=1:n %parfor бесмысленнен
            S_MEDD(2*ii-1) = I(ii);
            S_MEDD(2*ii) = I(ii);
            %tmp = [I(ii) I(ii)];
            %S_MEDD = [S_MEDD tmp];
        end
        %aend1=toc(ast1);
        %S_new = [5 4 5 9 3 8 9 7 6 2 7 1 2 8 1 6 10 3 10 4];
        %S_new_full = JeongFullS(S_new);
        %[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1);
        
        S0 = S_MEDD;
        %S=[];
        S = zeros(1,2*n);
        %S_new_par = [];
        %ri=1;
        %while ri<=2*n
        %ast1=tic;
        for  ri = 1:2*n   %этот цикл не параллелится!
            sr = S0(ri);
            sr2 = S0(ri+1:numel(S0));
            %T_min = inf;
            %S_min = [];
            %C_min = inf;
            T_S_new=ones(1,ri)*inf; %numel(S)
            SS1 = S(1:ri-1);
            parfor l=1:ri %parfor имеет смысл
                if l==1
                    S_new_par = [sr SS1];
                elseif l== ri
                    S_new_par = [SS1 sr];
                else
                    S_new_par = [SS1(1:l-1) sr SS1(l:ri-1)];
                end
                %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
                S_new = [S_new_par sr2];
                S_new_full = JeongFullS(S_new);
                %S_new_full = JeongFullSp(S_new);
                %важно! вместо d - d1
                %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
                % важно! добавил r
                [T_S_new(l),~,~,~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
            end
            [T_min, ind_m] = min(T_S_new);
            T_min_old = T_min;
            %S=S_min;
            if ind_m == 1
                %S = [sr S];
                S(2:ri) = S(1:ri-1);
                S(1) = sr;
            elseif ind_m == ri %numel(S) + 1
                %S = [S sr];
                S(ri) = sr;
            else
                %S = [S(1:ind_m-1) sr S(ind_m:numel(S))];
                S(ind_m+1:ri)=S(ind_m:ri-1);
                S(ind_m) = sr;
            end
            
            %конец MNEH начало FL
            if flag_1212_fl == 1
                rnd = rand();
                if rnd<Kross_fl
                    if ri>=3
                        %S_min = [];
                        T_S_new=ones(ri-1,ri)*inf; %numel(S)
                        SS1 = S(1:ri);
                        parfor ii=1:ri-1 %эффективен %numel(S)
                            %for ii=1:1
                            T_S_new(ii,:) = S_xchange(ii,SS1,ri,sr2, p1, n, m1, ma1, d1, su1, r);
                        end
                        [T_S_new1,ind1] = min(T_S_new);
                        [T_min,ind2] = min(T_S_new1);
                        i=ind1(ind2);
                        j=ind2;
                        if T_min<T_min_old
                            tempo = S(j);
                            S(j)=S(i);
                            S(i) = tempo;
                        end
                    end
                end
            end
        end
        %aend1=toc(ast1);
        if T_min<record
            xmin=S;
            record = T_min;
        end
        
        if record>0 && flag_branch ==1 %1111111
            U=xmin; %
            S_full=[];
            end_time_ex = [];
            %k2=1;
            %k3=1;
            %LocDom2=zeros(numel(U),2);
            %LocDom3=zeros(numel(U),3);
            % LocDom2 = [];
            % LocDom3 =[];
            
            Un=unique(U,'stable');
            %ast1 = tic; %
            %[LocDom2, LocDom3] = KUMZ_get_dom(S_full,Un, p,n,m,ma,d,su,end_time_ex);
            %[LocDom2, LocDom3] = KUMZ_get_domP(S_full,Un, p,n,m,ma,d,su,end_time_ex);
            [LocDom2, LocDom3] = KUMZ_get_domPsu(S_full,Un, p,n,m,ma,d,su,end_time_ex, r); %затратная по времени!
            %aend1 = toc(ast1);
            
            LocDom2t = unique(LocDom2,'rows');
            LocDom3t2 = unique(LocDom3,'rows');
            %LocDom3t2 = [];
            
            %U=1:n; %самый простой вариант
            %U=[U U]; %самый простой вариант
            S=[];
            %S_full=[];
            %end_time_ex = [];
            
            %условия доминирования 0 уровня то же пустые
            LocDom2 = [];
            LocDom3 = [];
            LocDom3t1 = [];
            
            ccount = 0;
            
            Cyc=2;
            %record1=record;
            
            %закроем пул ??? пока нет
            %matlabpool close
            if ~isempty(gcp('nocreate')) %нет
                delete(gcp('nocreate'));
            end
            
            %stt = tic;%
            %[~] = JeongBranchAndBound9(p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc);
            in_arg = {p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc,record,xmin,ccount};%,record,xmin
            %последовательный вариант
            %[~,xmin] = JeongBranchAndBound13(in_arg); % %важно! в данной функции сменилась модель хранения su
            
            %параллельный вариант
            [~,xmin] = JeongBranchAndBound11p(in_arg); % %важно! в данной функции сменилась модель хранения su
            
            %тест
            %[T_S_test,~,C2_S_new_test, end_t_test]= XieTargFun7r(p1, xmin, n, m1, ma1, d1, su1, r);
            
            %endt = toc(stt);
            S121_1 = xmin;
        else
            %S121_1 = [2 8 7 6 7 6 10 9 3 10 1 3 8 1 5 4 5 4 9 2];
            S121_1 = xmin;
        end
        
        %--------------------------------------------------------------------------
        % пробуем решить как хотели: (1-2-1-2)-2
        % в. 1. просто сложив (1-3-1)=1,(2,3,4)=2,(1-3-1)-2=1,(2,3,4)-2=2,
        % ma=[1 3 1 2 3 4 1 3 1 2 3 4 1 3 1 2 3 4 5];
        if flag_131 == 1
            % пробуем решить задачу (1,3,1)-1
            %         a1=p(7,:);
            %         b1=p(8,:);
            %         c1=p(9,:);
            % 7-8-9 уже решена смотри xmin2
            %
            %         a=sum(a1);
            %         b=sum(b1);
            %         c=sum(c1);
            %
            %         pj=[a1;b1]';
            %         S1= jhons(pj);
            %         Cab = WangTargFun(a1,b1,c1,S1);
            %         if Cab == a + c
            %             xmin1 = S1;
            %             %record1 = Cab;
            %             %exit=0;
            %         else
            %             F=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
            %             xmin1 = xmin;
            %             %record1 = record;
            %         end
            
            % пробуем решить задачу (1,3,1)-2
            a1=p(13,:);
            b1=p(14,:);
            c1=p(15,:);
            
            a=sum(a1);
            b=sum(b1);
            c=sum(c1);
            
            pj=[a1;b1]';
            S1= jhons(pj);
            Cab = WangTargFun(a1,b1,c1,S1);
            if Cab == a + c || flag_jhons ==1
                xmin3 = S1;
                %record2 = Cab;
                %exit=0;
            else
                [F,xmin3]=WangBranchAndBoundFirst(a1,b1,c1,S1,n, a, b, c);
                %xmin3 = xmin;
                %record2 = record;
            end
            % используя решения задач (1,3,1)-1 и (1,3,1)-2
            % модифицируем d
            [~,I]=sort(xmin2);
            d1=d(xmin2);
            
            for ii=1:n
                d1(ii)=d1(ii)-H/ii+H/n;
            end
            d2=d1(I);
            
            % модифицируем d по второй тройке
            [~,I]=sort(xmin3);
            d1=d2(xmin3);
            %H=200;
            
            for ii=1:n
                d1(ii)=d1(ii)-H/ii+H/n;
            end
            d2=d1(I);
        else
            d2=d;
        end
        
        p1=zeros(4,n);
        p1(1,:) = sum(p(7:9,:));
        p1(2,:) = sum(p(10:12,:));
        p1(3,:) = sum(p(13:15,:));
        p1(4,:) = sum(p(16:19,:));
        
        if sum(sum(p1))>0 % бывают случае когда второй части работ нет (1-2)
            %
            su1 = cell(1,n);
            %for ii=1:4 %parfor
            su1{1}(:,:)=su{7} + su{8} + su{9}; %sum(su{ii}(1:3,:));
            su1{2}(:,:)=su{10} + su{11} + su{12}; %sum(su{ii}(4:6,:));
            su1{3}(:,:)=su{13} + su{14} + su{15}; %sum(su{ii}(7:9,:));
            su1{4}(:,:)=su{16} + su{17} + su{18} + su{19}; %sum(su{ii}(10:12,:));
            %end
            
            %d1 = d - sum(p(1:6,:));
            %т.к. d уже с учетом реальных длительностей
            d1 = d2 - p(1,:) - p(2,:)*473 - p(3,:) - p(4,:)*18 - p(5,:)*473 - p(6,:); %d
            % важно для 2й задачи ограничение на r не уместно
            r1 = zeros(1,n);
            
            %применим алгоритм NEH для начального приближения
            S_NEH=1:n;
            m1=2;
            ma1=[1 2 1 2];
            %для NEH необх транспонировать
            %[F_NEH,NEH]=NawazHehCmax(p1',m1,ma1,S_NEH); %версия без r, т.к. вторая часть расписания от r не зависит
            [~,NEH]=NawazHehCmax2(p1',m1,ma1,S_NEH,r1); %
            
            %S_NEH = [];
            S_NEH = zeros(1,2*n);
            for ii=1:n %parfor бесмысленнен
                S_NEH(2*ii-1) = NEH(ii);
                S_NEH(2*ii) = NEH(ii);
                %tmp = [NEH(ii) NEH(ii)];
                %S_NEH = [S_NEH tmp]; %ускорение
            end
            
            S = JeongFullS(S_NEH);
            %S = JeongFullSp(S_NEH); %бесмысленнен
            %[T_NEH_1212,~,~,~]= XieTargFun7(p1, S, n, m1, ma1, d1, su1); %версия без r
            [T_NEH_1212,~,~,~]= XieTargFun7r(p1, S, n, m1, ma1, d1, su1, r1); %
            
            %%if T_NEH<record
            xmin=S_NEH;
            record = T_NEH_1212;
            clear T_NEH_1212;
            %%end
            
            %применим алгоритм модифицированный FL - MFL для начального приближения
            %MEDD
            d_m = d1 - (p1(3,:) + p1(4,:));
            [~, I] = sort(d_m);
            
            %S_MEDD = [];
            S_MEDD = zeros(1,2*n);
            for ii=1:n %parfor бесмысленнен
                S_MEDD(2*ii-1) = I(ii);
                S_MEDD(2*ii) = I(ii);
                %tmp = [I(ii) I(ii)];
                %S_MEDD = [S_MEDD tmp];
            end
            
            %S_new = [5 4 5 9 3 8 9 7 6 2 7 1 2 8 1 6 10 3 10 4];
            %S_new_full = JeongFullS(S_new);
            %[T_S_new,C1_S_new,C2_S_new, end_time_ex_S_new]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1);
            
            S0 = S_MEDD;
            %S=[];
            S = zeros(1,2*n);
            %S_new_par = [];
            %ri=1;
            %while ri<=2*n
            %ast1=tic;
            for  ri = 1:2*n   %этот цикл не параллелится!
                sr = S0(ri);
                sr2 = S0(ri+1:numel(S0));
                %T_min = inf;
                %S_min = [];
                %C_min = inf;
                T_S_new=ones(1,ri)*inf; %numel(S)
                SS1 = S(1:ri-1);
                for l=1:ri %parfor бесмысленнен
                    if l==1
                        S_new_par = [sr SS1];
                    elseif l== ri
                        S_new_par = [SS1 sr];
                    else
                        S_new_par = [SS1(1:l-1) sr SS1(l:ri-1)];
                    end
                    %здесь хорошо бы проверить условия доминирования, но алгоритм и так быстрый
                    S_new = [S_new_par sr2];
                    S_new_full = JeongFullS(S_new);
                    %S_new_full = JeongFullSp(S_new);
                    %важно! вместо d - d1
                    %[T_S_new,~,C2_S_new, ~]= XieTargFun7(p1, S_new_full, n, m1, ma1, d1, su1); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
                    % важно! добавил r
                    [T_S_new(l),~,~,~]= XieTargFun7r(p1, S_new_full, n, m1, ma1, d1, su1, r); % для 1-2-1-2 считается, что центры одинарные %(p1, S, n, m1, ma1, d, su1);
                end
                [T_min, ind_m] = min(T_S_new);
                T_min_old = T_min;
                %S=S_min;
                if ind_m == 1
                    %S = [sr S];
                    S(2:ri) = S(1:ri-1);
                    S(1) = sr;
                elseif ind_m == ri %numel(S) + 1
                    %S = [S sr];
                    S(ri) = sr;
                else
                    %S = [S(1:ind_m-1) sr S(ind_m:numel(S))];
                    S(ind_m+1:ri)=S(ind_m:ri-1);
                    S(ind_m) = sr;
                end
                
                %конец MNEH начало FL
                if flag_1212_fl == 1
                    rnd = rand();
                    if rnd<Kross_fl
                        if ri>=3
                            %S_min = [];
                            T_S_new=ones(ri-1,ri)*inf; %numel(S)
                            SS1 = S(1:ri);
                            parfor ii=1:ri-1 %эффективен numel(S)
                                %for ii=1:1
                                T_S_new(ii,:) = S_xchange(ii,SS1,ri,sr2, p1, n, m1, ma1, d1, su1, r);
                            end
                            [T_S_new1,ind1] = min(T_S_new);
                            [T_min,ind2] = min(T_S_new1);
                            i=ind1(ind2);
                            j=ind2;
                            if T_min<T_min_old
                                tempo = S(j);
                                S(j)=S(i);
                                S(i) = tempo;
                            end
                        end
                    end
                end
            end
            %aend1=toc(ast1);
            if T_min<record
                xmin=S;
                record = T_min;
            end
            clear T_min;
            clear T_min_old;
            
            if record>0 && flag_branch ==1 %22222222222
                U=xmin; %
                %S_full=[];
                end_time_ex = [];
                %k2=1;
                %$k3=1;
                %LocDom2=zeros(numel(U),2);
                %LocDom3=zeros(numel(U),3);
                
                %Un=unique(U,'stable');
                %[LocDom2, LocDom3] = KUMZ_get_dom(S_full,Un, p,n,m,ma,d,su,end_time_ex);
                %[LocDom2, LocDom3] = KUMZ_get_domP(S_full,Un, p,n,m,ma,d,su,end_time_ex);
                
                [LocDom2, LocDom3] = KUMZ_get_domPsu(S_full,Un, p,n,m,ma,d,su,end_time_ex); %затратная по времени!
                
                LocDom2t = unique(LocDom2,'rows');
                LocDom3t2 = unique(LocDom3,'rows');
                
                %U=1:n; %самый простой вариант
                %U=[U U]; %самый простой вариант
                S=[];
                %S_full=[];
                %end_time_ex = [];
                
                %условия доминирования 0 уровня то же пустые
                LocDom2 = [];
                LocDom3 = [];
                LocDom3t1 = [];
                
                Cyc=2;
                
                flag_branch=0;
                
                %закроем пулл
                if ~isempty(gcp('nocreate')) %нет
                    delete(gcp('nocreate'));
                end
                
                %stt1 = tic;%
                %[~] = JeongBranchAndBound9(p1,su1,r,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc);
                in_arg = {p1,su1,r1,d1,S,U,1,n,m1,ma1,end_time_ex,LocDom2,LocDom3,LocDom2t,LocDom3t1,LocDom3t2,Cyc,record,xmin,ccount};%,record,xmin
                %последовательный вариант
                %[~,xmin] = JeongBranchAndBound13(in_arg); % %важно! в данной функции сменилась модель хранения su
                
                %параллельный вариант
                [~,xmin] = JeongBranchAndBound11p(in_arg); %важно! в данной функции сменилась модель хранения su
                
                %endt1 = toc(stt1);
                S121_2 = xmin;
            else
                %S121_2 = [7 6 7 6 10 9 10 1 9 3 2 8 3 2 8 1 5 4 5 4];
                S121_2 = xmin;
            end
        else
            S121_2 = S121_1;
        end
        %--------------------далее базовые расписания S121_1-----------------------
        %пробуем простое решение
        S_1212f1 = KumzFullS12(S121_1,ma); % (perm) расписание - послед-й проход заказов:
        %S = KumzFullS12_2(S121_1); % (perm) расписание - послед-й по РЦ
        % важно! XieTargFun6 не параллелится
        %[T_1212f1,~,~, end_time_ex_1212f1]= XieTargFun7(p, S, n, m, ma, d, su);
        %[T_1212f1,~,~, end_time_ex_1212f1]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_1212f1,~,~, end_time_ex_1212f1]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_1212f1,~,~, end_time_ex_1212f1]= KUMZTargFun19f(p, S_1212f1, n, m, ma, d, su, mq, r, mq_u);
        FT_1212f1 = (T_1212f1 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212f1)) - Cmax_EDD);
        
        %для сравнения попробуем для 1-2 берется сортированный по a+b массив
        ab1 = sum(p(13:19,:));
        [~,S_ab] = sort(ab1);
        S_1212_ab=KumzFullS6_1(S121_1,S_ab);
        %[T_1212_ab,~,~, end_time_ex_1212_ab]= XieTargFun7(p, S, n, m, ma, d, su);
        %[T_1212_ab,~,~, end_time_ex_1212_ab]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_1212_ab,~,~, end_time_ex_1212_ab]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_1212_ab,~,~, end_time_ex_1212_ab]= KUMZTargFun19f(p, S_1212_ab, n, m, ma, d, su, mq, r, mq_u);
        FT_1212_ab = (T_1212_ab - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212_ab)) - Cmax_EDD);
        
        %--------------------далее базовые расписания S121_2-----------------------
        %пробуем простое решение
        S_1212f2 = KumzFullS12(S121_2,ma); % (perm) расписание - послед-й проход заказов:
        %S = KumzFullS12_2(S121_2); % (perm) расписание - послед-й по РЦ
        % важно! XieTargFun6 не параллелится
        %[T_1212f2,~,~, end_time_ex_1212f2] = XieTargFun7(p, S, n, m, ma, d, su);
        %[T_1212f2,~,~, end_time_ex_1212f2] = XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        %[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        %[T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_1212f2,~,~, end_time_ex_1212f2] = KUMZTargFun19f(p, S_1212f2, n, m, ma, d, su, mq, r, mq_u);
        FT_1212f2 = (T_1212f2 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212f2)) - Cmax_EDD);
        
        S_1212=KumzFullS6(S121_1,S121_2,ma);
        %[T_1212,~,~, end_time_ex_1212]= XieTargFun7(p, S, n, m, ma, d, su);
        %[T_1212,~,~, end_time_ex_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_1212,~,~, end_time_ex_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_1212,~,~, end_time_ex_1212]= KUMZTargFun19f(p, S_1212, n, m, ma, d, su, mq, r, mq_u);
        FT_1212 = (T_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_1212)) - Cmax_EDD);
        
        % для сравнения попробуем для 1-2 берется сортированный по a+b массив
        ab1 = sum(p(1:6,:));
        [~,S_ab] = sort(ab1);
        S_ab_1212=KumzFullS6(S_ab,S121_2,ma);
        %[T_ab_1212,~,~, end_time_ex_ab_1212]= XieTargFun7(p, S, n, m, ma, d, su);
        %[T_ab_1212,~,~, end_time_ex_ab_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_ab_1212,~,~, end_time_ex_ab_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_ab_1212,~,~, end_time_ex_ab_1212]= KUMZTargFun19f(p, S_ab_1212, n, m, ma, d, su, mq, r, mq_u);
        FT_ab_1212 = (T_ab_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_ab_1212)) - Cmax_EDD);
        
        % для сравнения попробуем первое 1-2 просто по Джонсону отсортировать
        a1=sum(p(1:3,:));
        b1=sum(p(4:6,:));
        p1=[a1;b1]';
        S_john= jhons(p1);
        S_joh_1212=KumzFullS6(S_john,S121_2,ma);
        %[T_joh_1212,~,~, end_time_ex_joh_1212]= XieTargFun7(p, S, n, m, ma, d, su);
        %[T_joh_1212,~,~, end_time_ex_joh_1212]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
        %[T_joh_1212,~,~, end_time_ex_joh_1212]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
        [T_joh_1212,~,~, end_time_ex_joh_1212]= KUMZTargFun19f(p, S_joh_1212, n, m, ma, d, su, mq, r, mq_u);
        FT_joh_1212 = (T_joh_1212 - T_EDD)*k1 + (1-k1)*(max(max(end_time_ex_joh_1212)) - Cmax_EDD);
        
        if flag_SACO2 == 1
            %-------------далее базовые расписания S121_1 - x - S121_2-----------------
            % в1. двойное перестановочное 2n - не сильно хорошо работает
            %[max_u_SACO,S_opt_SACO]=SACO_KUMZ_12_x_12(p, n, m, ma, d, su,S121_1,S121_2);  %для SACO нужно учесть r !!!
            % в2. одинарное перестановочное n(упрощение)
            %[max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_3(p, n, m, ma, d, su,S121_1,S121_2);  %
            % stt3 = tic;
            %не параллельная версия
            % [max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4(p, n, m, ma, d, su,S121_1,S121_2, r);  %
            % endt3 = toc(stt3);
            %stt3 = tic;
            %параллельная версия
            [~,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4p(p, n, m, ma, d, su,S121_1,S121_2, r,flag_fl3_2,flag_KUMZ, mq, mq_u, T_EDD, Cmax_EDD,k1);  %
            %endt3 = toc(stt3);
            %[max_u_SACO,S_opt_SACO,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_5(p, n, m, ma, d, su,S121_1,S121_2,mq);  %
            S_12_SACO_12=KumzFullS11(S1_s,S_opt_SACO,S2_e,S2_U);
            %[T_12_SACO_12,~,~, end_time_12_SACO_12]= XieTargFun7(p, S, n, m, ma, d, su);
            %[T_12_SACO_12,~,~, end_time_12_SACO_12]= XieTargFun12f(p, S, n, m, ma, d, su, mq);
            %[T_12_SACO_12,~,~, end_time_12_SACO_12]= KUMZTargFun17f(p, S, n, m, ma, d, su, mq, r, mq_u);
            [T_12_SACO_12,~,~, end_time_12_SACO_12]= KUMZTargFun19f(p, S_12_SACO_12, n, m, ma, d, su, mq, r, mq_u);
            FT_12_SACO_12 = (T_12_SACO_12 - T_EDD)*k1 + (1-k1)*(max(max(end_time_12_SACO_12)) - Cmax_EDD);
        else
            FT_12_SACO_12 = inf;
            S_12_SACO_12=S_joh_1212;
        end
        if flag_GA2 == 1
            %параллельная версия
             [~,S_opt_GA,S1_s,S2_e,S2_U ] = HGA_KUMZ( p, n, m, ma, d, su,S121_1,S121_2, r,mq, mq_u, flag_KUMZ,T_EDD, Cmax_EDD,k1);
            
            S_12_GA_12=KumzFullS11(S1_s,S_opt_GA,S2_e,S2_U);
            [T_12_GA_12,~,~, end_time_12_GA_12]= KUMZTargFun19f(p, S_12_GA_12, n, m, ma, d, su, mq, r, mq_u);
            FT_12_GA_12 = (T_12_GA_12 - T_EDD)*k1 + (1-k1)*(max(max(end_time_12_GA_12)) - Cmax_EDD);
        else
            FT_12_GA_12 = inf;
            S_12_GA_12=S_joh_1212;
        end
    else
        FT_1212_ab = inf;
        FT_1212f1 = inf;
        FT_1212f2 = inf;
        FT_12_SACO_12 = inf;
        FT_ab_1212 = inf;
        FT_joh_1212 = inf;
        FT_1212 = inf;
        S_1212_ab = S_EDD;
        S_1212f1 = S_EDD;
        S_1212f2 = S_EDD;
        S_12_SACO_12 = S_EDD;
        S_ab_1212 = S_EDD;
        S_joh_1212 = S_EDD;
        S_1212 = S_EDD;
    end
    %----------------------------------------------------------------------------------------------------
    %вывод лучшего результата в нужном формата
    FF = [FT_Maxp FT_NEH_f FT_1212_ab FT_1212f1 FT_1212f2 FT_12_SACO_12 FT_ab_1212 FT_joh_1212 FT_1212 FT_MNEH 0 FT_GA FT_SACO]; % 0 - EDD
    [~,ind] = min(FF);
    SS = [S_maxp; S_NEH_f; S_1212_ab; S_1212f1; S_1212f2; S_12_SACO_12; S_ab_1212; S_joh_1212; S_1212; S_MNEH; S_EDD; S_GA; S_SACO];
    Smin(glob_i,:) = SS(ind,:);
    
end

if n_glob>1
    
    S2min = [];
    for glob_i=1:n_glob
        Smin(glob_i,:) = Smin(glob_i,:) + n*(glob_i-1);
        S2min = [S2min Smin(glob_i,:)];
    end
    
    p=p_in(:,1:n*n_glob); %1:n
    
    %su=cell(1,19);
    for ii=1:19 %parfor % бессмысленнен
        su{ii}=zeros(n*n_glob,n*n_glob); %zeros(n_max,n_max)
    end
    
    % валки
    su{1}=su1_in(1:n*n_glob,1:n*n_glob);
    su{3}=su1_in(1:n*n_glob,1:n*n_glob);
    su{7}=su1_in(1:n*n_glob,1:n*n_glob);
    su{9}=su1_in(1:n*n_glob,1:n*n_glob);
    su{13}=su1_in(1:n*n_glob,1:n*n_glob);
    su{15}=su1_in(1:n*n_glob,1:n*n_glob);
    % ножи
    su{6}=su2_in{1}(1:n*n_glob,1:n*n_glob);
    su{12}=su2_in{2}(1:n*n_glob,1:n*n_glob);
    su{18}=su2_in{3}(1:n*n_glob,1:n*n_glob);
    
    %r=r_in(1:n*n_glob);
    %пока не уверен в правильности выгрузки r
    r = zeros(1,n*n_glob); %release dates - добавлены в алгоритм
    d = d_in(1:n*n_glob);
    [T_S2min,~,~, end_time_S2min]= KUMZTargFun19f(p, S2min, n*n_glob, m, ma, d, su, mq, r, mq_u);
    
else
    S2min = Smin;
end

%delete(gcp('nocreate'));

%-------доп оптимизация------------------------------
% версия с алгоритмом Xie
% в этом случае решаем задачу 1,3,1,(2,3,4),1 дважды
% для начала простая эвристика
%  2-я часть расписания строится по Wang 1-3-1
%  для 1-2 берется сортированный по a+b массив
% c1 = p(3,:);
% d1 = sum(p(4:6,:));
% e1 = p(7,:);
%
% c = sum(c1);
% d = sum(d1);
% e = sum(e1);
%
% pj=[c1;d1]';
% Sj= jhons(pj);
% F=WangBranchAndBoundFirst(c1,d1,e1,Sj,n, c, d, e);
%
% ab1 = a1 + b1;
% [ab1_sort,S_ab] = sort(ab1);

end