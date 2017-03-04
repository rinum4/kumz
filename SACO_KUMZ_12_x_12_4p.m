function [min_u_glob,xmin,S1_s,S2_e,S2_U]=SACO_KUMZ_12_x_12_4p(p, n, m, ma, d, su,S1_const,S2_const, r,flag_fl,flag_KUMZ, mq, mq_u, T_EDD, Cmax_EDD,k1)
%Ant Colony Optimization for
%ma=[(1 3 1)-S1 (2 3 4)-S1 (1 3 1)-x (2 3 4)-x (1 3 1)-S2 (2 3 4 5)-S2];
%важно! в данной функции сменилась модель хранения su

N_iter = 300; %максимальное количество итераций - классический алгоритм
N_max = N_iter/2; %максимальное количество итераций без изменения
M = ceil(n/2);   %количество муравьев ceil(n/3) 50

%Q = 1000000;          %параметр обучения 1
rho=0.9;       %параметр обучения 2 %коэфф испарения

%классическая реализация ACO
%nVar = n*numel(ma);
nVar=n;

%начинаем с S2 которое однозначно определено
[~,nc] = size(S2_const);
%S=[];
%S2_U=[];
S2_U=zeros(1,n);
S2_e=[]; %исправить нельзя т.к. формируется динамически
%упрощаем - x всегда вектор
%for i=nr:-1:1
cc=1;
if ~isempty(S2_const)
    for j=nc:-1:(nc-n+1) %parfor не пойдет т.к. S2_U
        ind = find(S2_U==S2_const(1,j), 1);
        if isempty(ind)
            cur_job = ones(1,7)*S2_const(1,j);
        else
            cur_job = ones(1,6)*S2_const(1,j);
        end
        %S2_U = [S2_U S2_const(1,j)];
        S2_U(cc) = S2_const(1,j);
        cc=cc+1;
        S2_e=[cur_job S2_e];
    end
end    
%end

%S1 - если дважды в S2_U пропуск
%S1_un=unique(S1,'stable');
[~,nc] = size(S1_const);
S1_s=[]; %исправить нельзя т.к. формируется динамически
%S1_U=[];
S1_U=zeros(1,n);
count=0;
cc=1;
%упрощаем - S всегда вектор
%for i=1:nr
if ~isempty(S1_const)
    for j=1:nc %parfor не пойдет т.к. S1_U
        if count<n
            ind = find(S2_U==S1_const(1,j));
            if numel(ind)==2
                ind2 = find(S1_U==S1_const(1,j),1);
                if isempty(ind2)
                    cur_job = ones(1,6)*S1_const(1,j);
                    %S1_U = [S1_U S1_const(1,j)];
                    S1_U(cc) = S1_const(1,j);
                    cc=cc+1;
                    S1_s=[S1_s cur_job];
                    count = count+1;
                else
                    %cur_job = [];
                end
            else
                cur_job = ones(1,6)*S1_const(1,j);
                %S1_U = [S1_U S1_const(1,j)];
                S1_U(cc) = S1_const(1,j);
                cc=cc+1;
                S1_s=[S1_s cur_job];
                count = count+1;
            end
        else
            break;
        end    
    end
end    
%end

if ~isempty(S2_const)
    T1 = unique(S2_const,'stable'); % работы которые необходимо спланировать
    %T1 = [T1 T1]; %
    
    [~,nc] = size(T1);
    %T2=[];
    T2 = zeros(1,nc);
    %for i=1:nr
    ne=1;
    for j=1:nc %parfor бесмысленнен
        ind1 = find(S1_U==T1(1,j));
        ind2 = find(S2_U==T1(1,j));
        if (numel(ind2)==2 && numel(ind1)==1) || (numel(ind2)==1 && numel(ind1)==2)
            %cur_job = [];
        elseif (numel(ind2)==2 && numel(ind1)==0) || (numel(ind2)==1 && numel(ind1)==1) || (numel(ind2)==0 && numel(ind1)==2)
            cur_job = T1(1,j);
            T2(ne) = cur_job;
            ne=ne+1;
        elseif (numel(ind2)==1 && numel(ind1)==0) || (numel(ind2)==0 && numel(ind1)==1)
            cur_job = [T1(1,j) T1(1,j)];
            T2(ne:ne+1) = cur_job;
            ne=ne+2;
        else
            cur_job = [T1(1,j) T1(1,j) T1(1,j)];
            T2(ne:ne+2) = cur_job;
            ne=ne+3;
        end
        %T2=[T2 cur_job];
        %T2() = cur_job;
    end
else
    %T2=1:n;
    [~, I] = sort(d);
    T2 = I;
end
%end

tau_max = 1/(1-rho); %начальная конентрация ферромонов 1*Q/(1-rho)
tau_min = tau_max/5;
tau=1*ones(nVar,nVar)*tau_max;   % матрица ферромонов 1*ones(nVar,nVar)

%S_opt = T2;
min_u_glob = inf;
xmin=T2;

iteration=1;
N_un = 0;
%fl_st=0;
S = zeros(M,n);
Cost = ones(1,M)*inf;
%for iteration=1:N % классический алгоритм
while N_un < N_max && iteration<N_iter % обрываем если функция стабилизировалась
    p0=log(iteration)/log(N_iter);
    
    %ast1=tic;
    parfor Ant = 1:M % parfor эффективен
        %не параллельная версия
        %[S(Ant,:), Cost(Ant)] = SACO_KUMZ_12_x_12_42p(tau,S1_s,T2,S2_e,S2_U,nVar,p0,p, n, m, ma, d, su, r);
        %параллельная версия
        [S(Ant,:), Cost(Ant)] = SACO_KUMZ_12_x_12_42pp(tau,S1_s,T2,S2_e,S2_U,nVar,p0,p, n, m, ma, d, su, r, flag_fl,flag_KUMZ, mq, mq_u, T_EDD, Cmax_EDD,k1);
    end %Ant
    %aend1=toc(ast1);
    
    [min_u,ind] = min(Cost(:));
    S_opt = S(ind,:); %оптимальное расписание среди муравьев
    
    %обновляем концентрацию ферромонов
    %в.1 T'Kindt
    % испарение
    tau=rho*tau;
    if iteration==1
        Q = min_u;
    end    
    for j=1:n
        %tau(S_opt(j),j) = tau(S_opt(j),j) + Q/min_u;
        tau(j,S_opt(j)) = tau(j,S_opt(j)) + Q/min_u;
    end
    
    %в.2.
    %     tour = zeros(1,n+1);
    %     for k=1:M
    %         tour(1:n)=S(k,:);
    %         %tour=[tour tour(1)]; %
    %         tour(n+1)=tour(1); %
    %         for l=1:nVar
    %             i=tour(l);
    %             j=tour(l+1);
    %             tau(i,j)=tau(i,j)+Q/Cost(k); %Q2
    %         end
    %     end
    %     % испарение
    %     tau=rho*tau;
    
    tau=max(tau,tau_min); % ограничиваем снижение ферромонов
    tau=min(tau,tau_max); % ограничиваем повышение ферромонов
    
    if min_u<min_u_glob % оптимум поменялся
        N_un = 0; 
        min_u_glob = min_u;
        xmin = T2(S(ind,:));
    else % оптимум не поменялся
        N_un = N_un+1; 
    end
    
    %для графика
    %if max_u_i(iteration)==0
    %    max_u_i(iteration) = max_u_i(iteration-1);
        %N_un = N_un+1; % оптимум не поменялся
    %end
    iteration = iteration + 1;
end %iteration
%min_u = 
%F=max_u;
%plot(max_u_i(:))
%xmin=xmin;
end %fcnt




