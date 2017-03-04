function [S, Cost] = SACO_KUMZ_12_x_12_42p(tau,S1_s,T2,S2_e,S2_U,nVar,p0,p, n, m, ma, d, su, r)
%L=T; от L в реализации ничего не зависит
S=zeros(1,nVar);
%S1=[];
%Cost(Ant)=0;

%if Ant>1 && fl_st==0
    tau_i=tau;
    position = 1;
    %max_ant=inf;
    %k_old=0;
    while position <= nVar
        pr = rand();
        if pr>=p0
            %диверсификация
            %выбираем случайно
            P=tau_i(:,position);
            P=P/sum(P); %приводим вероятности в интервал [0 1]
            
            rr=rand;
            C=cumsum(P);
            k=find(rr<=C,1,'first');
        else
            %интенсификация
            %выбираем по максимальной концентрации
            [~,k]=max(tau_i(:,position));
        end
        
        if isempty(k)
            %плохой муравей, начинаем сначала
            position = 1;
            tau_i=tau;
            %max_ant=inf;
            continue;
        end
        
        %S1(position)=k;
        
        S(position)=k;
        tau_i(k,:) = 0;
        %max_ant=LB;
        %L(k)=0; % от L в реализации ничего не зависит
        position = position + 1;
    end %pos
%else %хорошее нач приближение? - не в параллельной версии
%   S=1:nVar;
%   fl_st=1;
%end

%API - Adjacent Pairwise Interchange - попарная перестановка
%локальный поиск - перестановка 1-API,2-API,...,(n-j)-API

%т.к. в муравье расписание от 1 до nVar далее необходимо его
%перевести в норм вид
max_S = S;
full_S = KumzFullS11(S1_s,T2(S),S2_e,S2_U);
max_loc_s = XieTargFun7r(p, full_S, n, m, ma, d, su, r);

for j=1:nVar-1
    for i=j+1:nVar
        cur_S = S;
        tempo = cur_S(j);
        cur_S(j)=cur_S(i);
        cur_S(i) = tempo;
        full_S = KumzFullS11(S1_s,T2(cur_S),S2_e,S2_U);
        Cost_i_j = XieTargFun7r(p,full_S, n, m, ma, d, su, r);
        if Cost_i_j<max_loc_s
            max_S = cur_S;
            max_loc_s = Cost_i_j;
        end
    end
end

%вписываем в муравья лучшее локальное перестановочное расписание
S =  max_S;
%Cost(Ant) = CroceTargFun(p,S);
Cost = max_loc_s;

%if Cost < max_u %в параллельной версии max_u не может передаваться между
%итерациями !
%     max_u=Cost;
%     %max_u_i(iteration)=Cost(Ant); %для графика
%     xmin = T(S);
%     %N_un = 0; % оптимум поменялся
%     %else
%     %    N_un = N_un+1; % оптимум не поменялся
%end

end

