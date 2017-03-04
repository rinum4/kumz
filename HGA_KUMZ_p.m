function  [Child, Cost] = HGA_KUMZ_p(Par1,Par2,Kross_p,Mut_p,Other_p, S1_s,T2,S2_e,S2_U,p, n, m, ma, d, su, r, mq, mq_u, flag_KUMZ,T_EDD, Cmax_EDD,k1)
%функция генерации поколений и расчета целевой функции

Child = Par1;

%кроссовер двухточечный 
%Murata(1996)
rnd = rand();
if rnd<Kross_p
    pos1 = ceil(rand()*n);
    pos2 = ceil(rand()*n);
    
    if pos2<pos1
        tmp=pos1;
        pos1=pos2;
        pos2=tmp;
    elseif pos1==pos2
        if pos2==n
            pos1 = pos2 - 1;
        else
            pos2 = pos1+1;
        end
    end    
    
    Child(pos1:pos2) = Par1(pos1:pos2); % выбираем середину из 1го родителя
    pos=1;
    %выбираем остальные элементы из второго родителя
    for i=1:pos1
        for j=pos:n
            ind = find(Child == Par2(j),1);
            if ~isempty(ind)
                pos=pos+1;
                continue;
            else
                Child(i) = Par2(j);
                break;
            end
        end
    end
    for i=pos2:n
        for j=pos:n
            ind = find(Child == Par2(j),1);
            if ~isempty(ind)
                pos=pos+1;
                continue;
            else
                Child(i) = Par2(j);
                break;
            end
        end
    end
end

%мутация сдвигом
%Murata(1996)
rnd = rand();
if rnd<Mut_p
    pos1 = ceil(rand()*n);
    pos2 = ceil(rand()*n);
    if pos2<pos1
        tmp=pos1;
        pos1=pos2;
        pos2=tmp;
    elseif pos1==pos2
        if pos2==n
            pos1 = pos2 - 1;
        else
            pos2 = pos1+1;
        end
    end    
    
    tmp = Child(pos2);
    for i = pos2:-1:pos1+1
        Child(i)=Child(i-1);
    end
    Child(pos1)=tmp;
end    

%гибридный поиск подзадачи с использованием NEH
rnd = rand();
if rnd<Other_p
    pos1 = ceil(rand()*n);
    pos2 = ceil(rand()*n);
    if pos2<pos1
        tmp=pos1;
        pos1=pos2;
        pos2=tmp;
    elseif pos1==pos2
        if pos2==n
            pos1 = pos2 - 1;
        else
            pos2 = pos1+1;
        end
    end
    
    np = pos2-pos1+1;
    if np>1
        p1=zeros(numel(ma),np);
        %NEH считается для подзадачи с pos1 до pos2
        for i=1:np
            p1(:,i)=p(:,Child(i+pos1-1));
        end
        S=1:np;
        [~,NEH]=NawazHehCmax2(p1',m,ma,S,r);
        %[~,NEH]=NawazHehCmax3(p1',m,ma,S,r);
        
        tmp = Child(pos1:pos2);
        Child(pos1:pos2) = tmp(NEH);
    end
end

if ~isempty(S2_e)
    full_S = KumzFullS11(S1_s,T2(Child),S2_e,S2_U);
else
    full_S = KumzFullS(T2(Child),ma);
end    
if flag_KUMZ == 0
    %быстрый вариант без учета параллельных РЦ
    [Cost,~,Cmax2,~] = XieTargFun7r(p, full_S, n, m, ma, d, su, r);
    Cost = (Cost - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
else
    %полный вариант
    [Cost,~,Cmax2,~] = KUMZTargFun19f(p, full_S, n, m, ma, d, su, mq, r, mq_u); 
    Cost = (Cost - T_EDD)*k1 + (1-k1)*(Cmax2 - Cmax_EDD);
end   
end

