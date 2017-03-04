function [SC]=NawazTargFun2(p,m,ma,S,r)
%важно! в данной функции начинаем использовать r

SC = 0;
n_s = numel(S);
n_ma = numel(ma);
p1 = p(S,:);
r1 = r(S);
h=zeros(n_s, n_ma); %начало iй работы на j-м РЦ
e=zeros(m,1); %конец  послед операции на РЦ

%расписание перестановочное
%так же считаем, что количество машин равно n_ma (развертываем для общности)
for j=1:n_ma
    if j==1
        for i=1:n_s
            if i==1
                %h(i,j)= 0;
                h(i,j)= r1(i);
            else
                h(i,j) = max(h(i-1,j),r1(i))+p1(i-1,j);
            end
        end %i
        e(1)=h(i,j) + p1(i,j);
    else
        k=ma(j);
        
        for i=1:n_s
            if i==1
                h(i,j)= max(h(i,j-1) + p1(i,j-1), e(k));
            else
                h(i,j)= max(h(i-1,j)+p1(i-1,j), h(i,j-1)+p1(i,j-1));
            end
        end % i
        e(k)=h(i,j) + p1(i,j);
    end %j==1
end %j
SC = h(i,j) + p1(i,j);
end