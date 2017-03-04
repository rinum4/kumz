function [F_max, a_max]=NawazHehCmax3(p,m,ma,S,r)
%важно! в данной функции начинаем использовать r
%важно! версия функции с полным расчетом
% KUMZTargFun19f вместо NawazTargFun2

%global F_max;
%global a_max;

p1 = p(S,:);
n1=numel(S);
r1 = r(S);

%n_ma = numel(ma);

%в начале сортируем по возрастанию времени выполнения
%t = sum(p1');
t = sum(p1,2);
[~, I]=sort(t);

%p_s = p(I,:);
a(1)=I(1);
a(2)=I(2);
F1 = KUMZTargFun19f(p, [a(1) a(2)], 2, m, ma, d, su, mq, r, mq_u);
F2 = KUMZTargFun19f(p, [a(2) a(1)], 2, m, ma, d, su, mq, r, mq_u);

[F, I2] = min([F1 F2]);
%[F, I2] = min([NawazTargFun2(p1,m,ma,[a(1) a(2)], r1) NawazTargFun2(p1,m,ma,[a(2) a(1)], r1)]);

if I2~=1
    a(1)=I(2);
    a(2)=I(1);
end    

if n1>2
    F_max=inf;
    
    for i=3:n1
        for j=1:numel(a)+1
            if j==1
                a1=[I(i) a];
                F_max = %NawazTargFun2(p1,m,ma,a1, r1);
                a_max = a1;
            else
                a1=[a(1:j-1) I(i) a(j:numel(a))];
                F = NawazTargFun2(p1,m,ma,a1, r1);
                if F<F_max
                    F_max=F;
                    a_max = a1;
                end
            end   %j==1
        end  %j
        a=a_max;
    end
else
    F_max=F;
    a_max=a;
end
%вернем в исходных индексах
a_max = S(a_max);
end