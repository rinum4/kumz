function [SC]=NawazTargFun(p,m,ma,S)

SC = 0;
n_s = numel(S);
n_ma = numel(ma);
p1 = p(S,:);
h=zeros(n_s, n_ma); %начало iй работы на j-м РЦ
e=zeros(m,1); %конец  послед операции на РЦ

%расписание перестановочное
%так же считаем, что количество машин равно n_ma (развертываем для общности)
for j=1:n_ma
    if j==1
        for i=1:n_s
            if i==1
                h(i,j)= 0;
            else
                h(i,j)= h(i-1,j)+p1(i-1,j);
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