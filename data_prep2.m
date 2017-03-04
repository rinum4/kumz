function [id,id_k] = data_prep2(db_data, id, nr) %db_data
%ориентируемся на заказы

% global db_data;
% global id;
% global nr;

% p=0;
% su=0;

max_id=inf; %30

id_k=1;

for i=1:nr
    if id(i)==0 && ~isnan(db_data{i,10})
        id(i) = id_k;
        
        if isnan(db_data{i,31})
            break;
        end
        
        id =   search_val2(db_data{i,31},i+1,id_k, db_data, id, nr);   
        
        id_k=id_k + 1;
        %         if id_k==51
        %             stop=1;
        %         end
        if id_k>max_id
            break;
        end    
    end    
end

end

function [id] = search_val2(val, i_s, id_k, db_data, id, nr)

% global db_data;
% global id;
% global nr;

max_i = inf; %5000 inf

%tmp=0;

for i=i_s:nr
    if i-i_s > max_i
        break;
    end    
    if strcmp(db_data{i,31}, val) && id(i)==0 && ~isnan(db_data{i,10})
        id(i) = id_k;
        if isnan(db_data{i,31})
            break;
        end    
        id =   search_val2(db_data{i,31}, i, id_k, db_data, id, nr);   
    end    
end    

end