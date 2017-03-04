function [ Costs ] =  SACO_KUMZ_12_x_12_42ppp(k,S,S1_s,T2,S2_e,S2_U,nVar,p, n, m, ma, d, su, r, flag_KUMZ, mq, mq_u,T_EDD, Cmax_EDD,k1)
Costs=ones(1,nVar)*inf;
for j=k+1:nVar %parfor %тест показал, что внутр пар цикл не эффективен
    cur_S = S;
    tempo = cur_S(k);
    cur_S(k)=cur_S(j);
    cur_S(j) = tempo;
    if ~isempty(S2_e)
        full_S = KumzFullS11(S1_s,T2(cur_S),S2_e,S2_U);
    else
        full_S = KumzFullS(T2(cur_S),ma);
    end    
    if flag_KUMZ == 0
        [min_loc_s,~,Cmax2,~] = XieTargFun7r(p, full_S, n, m, ma, d, su, r);
        %Costs(j) = XieTargFun7r(p,full_S, n, m, ma, d, su, r);
        Costs(j) = (min_loc_s - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
    else
        [min_loc_s,~,Cmax2,~] = KUMZTargFun19f(p, full_S, n, m, ma, d, su, mq, r, mq_u); 
        %Costs(j) = KUMZTargFun19f(p, full_S, n, m, ma, d, su, mq, r, mq_u); 
        Costs(j) = (min_loc_s - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
    end   
end

end

