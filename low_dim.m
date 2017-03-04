function [id_t_o2,id_in_o2,r_in_o2,d_in_o2,p_in_o2,su1_in_o2,su2_in_o2,S_out,S_umerge,S_merge] = low_dim(id_t,id_in,r_in,d_in,p_in,su1_in,su2_in)
% снизим размерность задачи!

%1. Если в работе только 1 задачи и она склад(2) или финишная операция(19) - работу выкидываем из
%рассмотрения, сохраняя массив соотв-ий. перед выводом все такие задачи
%ставим в начало...
%2. Если в работе только 2 задачи и это склад и финиш -  работу выкидываем из
%рассмотрения

[nr,nc]=size(p_in);
p_in_o = zeros(nr,nc); %+
id_t_o = zeros(1,nc);   %+
id_in_o = zeros(1,nc); %+
r_in_o = zeros(1,nc);   %+
d_in_o = zeros(1,nc);  %+
% su1_in_o = zeros(nc,nc);
% su2_in_o = cell(1,3);
% su2_in_o{1}=zeros(nc,nc);
% su2_in_o{2}=zeros(nc,nc);
% su2_in_o{3}=zeros(nc,nc);
su1_in_o = su1_in;
su2_in_o = su2_in;

S_out = zeros(1,nc);
S_in = zeros(1,nc);

k1=1;
k2=1;
for i=1:nc %цикл по заказам
    %     if i==212
    %         stop=1;
    %     end
    n=find(p_in(:,i));
    if numel(n)==1 && (n==2 || n==19) 
        S_out(k2)=i;
        
        %важно: сместить элементы матриц проще, чем циклами по S_out вставлять
        %значения в новую матрицу!
        su1_in_o(:,i:nc - k2-1) = su1_in_o(:,i+1:nc - k2);
        su1_in_o(i:nc - k2-1,:) = su1_in_o(i+1:nc - k2,:); %su - матрица квадратная
        
        su2_in_o{1}(:,i:nc - k2-1) = su2_in_o{1}(:,i+1:nc - k2);
        su2_in_o{1}(i:nc - k2-1,:) = su2_in_o{1}(i+1:nc - k2,:); 
        su2_in_o{2}(:,i:nc - k2-1) = su2_in_o{2}(:,i+1:nc - k2);
        su2_in_o{2}(i:nc - k2-1,:) = su2_in_o{2}(i+1:nc - k2,:); 
        su2_in_o{3}(:,i:nc - k2-1) = su2_in_o{3}(:,i+1:nc - k2);
        su2_in_o{3}(i:nc - k2-1,:) = su2_in_o{3}(i+1:nc - k2,:); 
        
        k2 = k2 + 1;
    elseif numel(n)==2 && (n(1)==2 && n(2)==19) 
        S_out(k2)=i;
        
        su1_in_o(:,i:nc - k2-1) = su1_in_o(:,i+1:nc - k2);
        su1_in_o(i:nc - k2-1,:) = su1_in_o(i+1:nc - k2,:); %su - матрица квадратная
        
        su2_in_o{1}(:,i:nc - k2-1) = su2_in_o{1}(:,i+1:nc - k2);
        su2_in_o{1}(i:nc - k2-1,:) = su2_in_o{1}(i+1:nc - k2,:);
        su2_in_o{2}(:,i:nc - k2-1) = su2_in_o{2}(:,i+1:nc - k2);
        su2_in_o{2}(i:nc - k2-1,:) = su2_in_o{2}(i+1:nc - k2,:);
        su2_in_o{3}(:,i:nc - k2-1) = su2_in_o{3}(:,i+1:nc - k2);
        su2_in_o{3}(i:nc - k2-1,:) = su2_in_o{3}(i+1:nc - k2,:);
        
        k2 = k2 + 1;
    else
        S_in(k1)=i;
        p_in_o(:,k1) = p_in(:,i);
        id_t_o(k1) = id_t(i);
        id_in_o(k1) = id_in(i);
        r_in_o(k1) = r_in(i);
        d_in_o(k1) = d_in(i);
        
                
        k1=k1+1;
    end    
end    
%снизим размерность
%[i,j] = find(p_in_o);
S_in = S_in(1:k1-1);
p_in_o = p_in_o(:,1:k1-1);
%p_in_o = p_in_o(1:max(i),:);
id_t_o = id_t_o(1:k1-1);
id_in_o = id_in_o(1:k1-1);
r_in_o = r_in_o(1:k1-1);
d_in_o = d_in_o(1:k1-1);

su1_in_o = su1_in_o(:,1:k1-1);
su1_in_o = su1_in_o(1:k1-1,:);

su2_in_o{1} = su2_in_o{1}(:,1:k1-1);
su2_in_o{1} = su2_in_o{1}(1:k1-1,:);
su2_in_o{2} = su2_in_o{2}(:,1:k1-1);
su2_in_o{2} = su2_in_o{2}(1:k1-1,:);
su2_in_o{3} = su2_in_o{3}(:,1:k1-1);
su2_in_o{3} = su2_in_o{3}(1:k1-1,:);

%[i,j] = find(S_out);
S_out = S_out(1:k2-1);

%3. Если работы совпадают, то они планируются как одна
% - совпадает date_req +/- 2%
% - совпадают операции и их длительность +/- 2%
% - совпадает оснастка
[nr,nc]=size(p_in_o);

p_in_o2 = zeros(nr,nc); %+
id_t_o2 = zeros(1,nc);   %+
id_in_o2 = zeros(1,nc); %+
r_in_o2 = zeros(1,nc);   %+
d_in_o2 = zeros(1,nc);  %+

su1_in_o2 = su1_in_o; %+
su2_in_o2 = su2_in_o; %+

S_u = zeros(1,nc);
S_umerge=zeros(1,nc);
S_merge=zeros(2,nc);
k1=1;
k2=1;
for i=1:nc %цикл по оставшимся заказам
    %     if i==164
    %         stop=1;
    %     end
    ind = find(S_u==i,1);
    if ~isempty(ind)
        continue;
    end    
   
    S_umerge(k1)=S_in(i);
    p_in_o2(:,k1) = p_in_o(:,i);
    id_t_o2(k1) = id_t_o(i);
    id_in_o2(k1) = id_in_o(i);
    r_in_o2(k1) = r_in_o(i);
    d_in_o2(k1) = d_in_o(i);
       
    for j=i+1:nc %цикл по другим заказам
        if abs(d_in_o(i)/d_in_o(j)-1)<0.02 || (d_in_o(j)==0 && d_in_o(i)==0)
            if abs(sum(p_in_o(:,i))/sum(p_in_o(:,j))-1)<0.02
                op1=find(p_in_o(:,i));
                op2=find(p_in_o(:,j));
                if numel(op1) == numel(op2)
                    if sum(op1==op2)==numel(op1) %операции совпадают
                        if su1_in_o(i,j)==0 % валки совпадают
                            if su2_in_o{1}(i,j)==0 && su2_in_o{2}(i,j)==0 && su2_in_o{3}(i,j)==0  %  ножи совпадают
                                S_merge(1,k2)=S_in(i);
                                S_merge(2,k2)=S_in(j);
                                S_u(k2)=j;
                                
                                %добавим операционное время
                                p_in_o2(:,k1) = p_in_o2(:,k1) + p_in_o(:,j);
                                
                                %важно: сместить элементы матриц проще, чем циклами по S_merge вставлять
                                %значения в новую матрицу!
                                su1_in_o2(:,j:nc - k2-1) = su1_in_o2(:,j+1:nc - k2);
                                su1_in_o2(j:nc - k2-1,:) = su1_in_o2(j+1:nc - k2,:); %su - матрица квадратная
                                
                                su2_in_o2{1}(:,j:nc - k2-1) = su2_in_o2{1}(:,j+1:nc - k2);
                                su2_in_o2{1}(j:nc - k2-1,:) = su2_in_o2{1}(j+1:nc - k2,:);
                                su2_in_o2{2}(:,j:nc - k2-1) = su2_in_o2{2}(:,j+1:nc - k2);
                                su2_in_o2{2}(j:nc - k2-1,:) = su2_in_o2{2}(j+1:nc - k2,:);
                                su2_in_o2{3}(:,j:nc - k2-1) = su2_in_o2{3}(:,j+1:nc - k2);
                                su2_in_o2{3}(j:nc - k2-1,:) = su2_in_o2{3}(j+1:nc - k2,:);
                                
                                k2=k2+1;
                            end
                        end
                    end
                end
            end    
        end    
    end  
    
    k1=k1+1;
end    
%снизим размерность
[~,j] = find(S_merge);
S_merge = S_merge(:,1:max(j));
S_umerge = S_umerge(1:k1-1);
p_in_o2 = p_in_o2(:,1:k1-1);
%p_in_o = p_in_o(1:max(i),:);
id_t_o2 = id_t_o2(1:k1-1);
id_in_o2 = id_in_o2(1:k1-1);
r_in_o2 = r_in_o2(1:k1-1);
d_in_o2 = d_in_o2(1:k1-1);

su1_in_o2 = su1_in_o2(:,1:k1-1);
su1_in_o2 = su1_in_o2(1:k1-1,:);

su2_in_o2{1} = su2_in_o2{1}(:,1:k1-1);
su2_in_o2{1} = su2_in_o2{1}(1:k1-1,:);
su2_in_o2{2} = su2_in_o2{2}(:,1:k1-1);
su2_in_o2{2} = su2_in_o2{2}(1:k1-1,:);
su2_in_o2{3} = su2_in_o2{3}(:,1:k1-1);
su2_in_o2{3} = su2_in_o2{3}(1:k1-1,:);
end

