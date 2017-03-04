function  [db_data1,id_t,id_in,r_in,d_in,p_in,su1_in,su2_in] = test_data_prep5(db_data)
%db_data_a = table2array(db_data);

%db_data = sortrows(db_data,'link_ingr','ascend');
flag_id = 1;
flag_r = 1;
flag_d = 1;
flag_p = 1;
flag_su1 = 1;
flag_su2 = 1;

% global db_data;
% global id;
% global nr;

%db_data = db_data_in;
%n=120;
%
[nr,~]=size(db_data);
%

if flag_id==1
    id_in=zeros(nr,1);

    db_data = sortrows(db_data,'link_ingr','ascend');%sortrows(db_data,'link_ingr','ascend');
    db_data1 = table2cell(db_data(:,1:31));
    [id_in,id_k] = data_prep(db_data1, id_in, nr);
    %[id_in,id_k] = data_prep2(db_data1, id_in, nr);
    [db_data, ind] = sortrows(db_data,'date_req','ascend');%sortrows(db_data,'date_req','ascend');
    % id = data_prep2(db_data, id, nr);
    id_t = id_in;
    id_in = id_in(ind);
    db_data = sortrows(db_data,'link_ingr','ascend');%sortrows(db_data,'link_ingr','ascend');
    % id = restore(db_data, id, nr);
end

n=id_k-1;
%n=10;

st_t = now;
%date_b = datenum(1970,01,01);
%date_b = datenum(1971,01,01); %+1год
date_b = datenum(1970,01,01); % +6мес
id_in_u = unique(id_in,'stable');
id_in_u = nonzeros(id_in_u);%find(id_in_u~=0);

if flag_r==1
    r_in=zeros(1,n);
    
    for i=1:n %цикл по заказам
        f_srch = id_in_u(i);
        if f_srch == 0
            continue;
        end
        id_n = find(id_t==f_srch);
        %stage=0;
        for j=1:numel(id_n) %цикл по строкам
            jj= id_n(j);
            %if id_in(jj)==f_srch % i
            if ~isnan(db_data1{jj,12})
                date_r = date_b + db_data1{jj,12}/3600/24;
                r_in(i) = max((date_r - st_t)*24*60,0);
                break;
            end
            %end
        end
    end
end

if flag_d==1
    d_in=zeros(1,n);
    
    for i=1:n %цикл по заказам
        f_srch = id_in_u(i);
        if f_srch == 0
            continue;
        end
        id_n = find(id_t==f_srch);
        %stage=0;
        for j=1:numel(id_n) %цикл по строкам
            jj= id_n(j);
            %if id_in(jj)==f_srch % i
            if ~isnan(db_data1{jj,22})
                date_r = date_b + db_data1{jj,22}/3600/24;
                d_in(i) = (date_r - st_t)*24*60;
                break;
            end
            %end
        end
    end
end

%---------------------------------------------------------
if flag_p==1
    p_in = zeros(19,n);
    for i=1:n %цикл по заказам
        f_srch = id_in_u(i);
        if f_srch == 0
            continue;
        end
        if f_srch == 151
            stop=1;
        end
        stage=0;
        fl1=0;
        fl2=0;
        fl3=0;
        id_n = find(id_t==f_srch);
        for j=1:numel(id_n) %цикл по строкам
            %if id_n(j)==f_srch % i
            %игнорим транспортные операции и операции хранения
            %важно! игнрить эти операции нельзя
            %if db_data1{j,4}{1}=='K' || db_data1{j,4}{1}=='T'
            %    continue;
            %end
            jj=id_n(j);
            
            %умножать на количество как оказалось не надо!
            %if ~isnan(db_data1{jj,20})
            %    durat = db_data1{jj,16}/60*db_data1{jj,20};
            %else
            durat = db_data1{jj,16}/60;
            %end
            
            if strcmp(db_data1{jj,3},'S01::ST')
                %p(1,i)=p(1,i) + db_data1{j,16}/60;
                continue;
            end % S01::ST
            
            if strcmp(db_data1{jj,3},'CRM::RS') ||  strcmp(db_data1{jj,3},'CRM::TX') %1рц
                if stage == 1 || stage == 5 || stage == 11
                    %if db_data1{j,4}~='K' && db_data1{j,4}~='T'
                    stage = stage + 1;
                    %end
                end
                
                if strcmp(db_data1{jj,3},'CRM::TX')
                    durat = durat/14; %14 тележек
                end
                
                switch stage
                    case 0
                        p_in(1,i)=p_in(1,i) + durat;
                    case 2
                        p_in(3,i)=p_in(3,i) +durat;
                        fl1=1;
                    case 6
                        p_in(7,i)=p_in(7,i) + durat;
                    case 8
                        p_in(9,i)=p_in(9,i) + durat;
                        fl2=1;
                    case 12
                        p_in(13,i)=p_in(13,i) + durat;
                    case 14
                        p_in(15,i)=p_in(15,i) + durat;
                        fl3=1;
                end
                continue;
            end    %'CRM::RS'
            
            if strcmp(db_data1{jj,3},'CTI::ST')  ||  strcmp(db_data1{jj,3},'CTI::TX') %473 ячейки
                if (fl1==1 && stage==2) || (fl2==1 && stage==8) || (fl3==1 && stage==14)
                  continue; %выход склад после 2го прохода не обязательный
                end
                if stage == 0 || stage == 6 || stage == 12
                    %if db_data1{j,4}~='K' && db_data1{j,4}~='T'
                    stage = stage + 1;
                    %end
                end
                
                switch stage
                    %case 0
                    %p_in(1,i)=p_in(1,i) + durat; %относится к 1-й операции
                    case {1,2}
                        p_in(2,i)=p_in(2,i) + durat/473; %10800 144000
                    case {4,5}
                        p_in(5,i)=p_in(5,i) + durat/473; %36000 144000
                    case {7,8}
                        p_in(8,i)=p_in(8,i) + durat/473; %10800 144000
                    case {10,11}
                        p_in(11,i)=p_in(11,i) + durat/473; %36000 144000
                    case {13,14}
                        p_in(14,i)=p_in(14,i) + durat/473; %10800 144000
                    case {16,17}
                        p_in(17,i)=p_in(17,i) + durat/473; %36000 144000
                end
                if stage~=2 && stage~=8 && stage~=14
                    if db_data1{jj,4}~='K' && db_data1{jj,4}~='T'
                        if stage==1 || stage==4 || stage==7 || stage==10 || stage==13 || stage==16 %если строки повторяются - не увеличивать
                            stage = stage + 1;
                        end
                    end
                end
                continue;
            end %'CTI::ST'
            
            if strcmp(db_data1{jj,3},'EBC::FX')  ||  strcmp(db_data1{jj,3},'EBC::TX') %18 рц
                if stage == 0 || stage == 6 || stage == 12
                    %if db_data1{j,4}~='K' && db_data1{j,4}~='T'
                    stage = stage + 3;
                    %end
                end
                
                if stage == 2 || stage == 8 || stage == 14
                    %if db_data1{j,4}~='K' && db_data1{j,4}~='T'
                    stage = stage + 1;
                    %end
                end
                
                if stage<3
                    stage=3;
                end
                
                switch stage
                    case {3,4}
                        p_in(4,i)=p_in(4,i) + durat/18;
                    case {9,10}
                        p_in(10,i)=p_in(10,i) +  durat/18;
                    case {15,16}
                        p_in(16,i)=p_in(16,i) +  durat/18;
                end
                if db_data1{jj,4}~='K' && db_data1{jj,4}~='T'
                    if stage==3 || stage==10 || stage==15 %если строки повторяются - не увеличивать
                        stage = stage + 1;
                    end
                end
                continue;
            end    %'EBC::FX'
            
            if strcmp(db_data1{jj,3},'G03::CU')  ||  strcmp(db_data1{jj,3},'G03::IN')  ||  strcmp(db_data1{jj,3},'G03::PA')  ||  strcmp(db_data1{jj,3},'G03::PK')  %1рц
                if stage<5
                    stage=5;
                end
                switch stage
                    case {5,6} %7 под вопросом
                        p_in(6,i)=p_in(6,i) + durat;
                    case {11,12}
                        p_in(12,i)=p_in(12,i) + durat;
                    case {17,18}
                        p_in(18,i)=p_in(18,i) + durat;
                end
                if db_data1{jj,4}~='K' && db_data1{jj,4}~='T'
                    if stage==5 || stage==11 || stage==17 %если строки повторяются - не увеличивать
                        stage = stage + 1;
                    end
                end
                continue;
            end    %'G03::CU'
            %2 рц
            if strcmp(db_data1{jj,3},'G01::CU')  ||  strcmp(db_data1{jj,3},'G01::IN')  ||  strcmp(db_data1{jj,3},'G01::PA')  ||  strcmp(db_data1{jj,3},'G01::PK') || ...
                    strcmp(db_data1{jj,3},'G02::CU')  ||  strcmp(db_data1{jj,3},'G02::IN')  ||  strcmp(db_data1{jj,3},'G02::PA')  ||  strcmp(db_data1{jj,3},'G02::PK')
                if stage<19
                    stage=19;
                end
                p_in(19,i)=p_in(19,i) + durat/2;
                if db_data1{jj,4}~='K' && db_data1{jj,4}~='T'
                    if stage<19 %если строки повторяются - не увеличивать
                        stage = stage + 1;
                    end
                end
                continue;
            end
            %end
        end %j цикл по строкам
        %нужна ли эта заплатка??
        %     if stage <18
        %         p(19,i)=p(19,i) + 40/2;
        %     end
    end
end
%----------------------------------------------------------------------
if flag_su1==1
    su1_in=zeros(n,n);
    
    for i=1:n %цикл по заказам
        f_srch = id_in_u(i);
        if f_srch == 0
            continue;
        end
        %      if f_srch == 60
        %         stop=1;
        %      end
        
        id_n = find(id_t==f_srch);
        
        search_val = 'NaN';
        %cyl=-1;
        for l=1:numel(id_n) %цикл по строкам
            ll=id_n(l);
            %if id_in(ll)==f_srch % i
            if strcmp(db_data1{ll,3},'CRM::RS')
                search_val = db_data1{ll,28};
                break;
            end
            %end
        end %l
        if ~strcmp(search_val, 'NaN')
            for jj=i+1:n %цикл по заказам
                f_srch = id_in_u(jj);
                id_n = find(id_t==f_srch);
                
                for l=1:numel(id_n) %цикл по строкам
                    ll=id_n(l);
                    %if id_in(ll)==f_srch % j
                    if strcmp(db_data1{ll,3},'CRM::RS')
                        if ~strcmp(db_data1{ll,28},search_val)
                            su1_in(i,jj)=60;
                            su1_in(jj,i)=60;
                        end
                        break;
                    end
                    %end
                end % l
            end % jj
        end
    end % i
end

if flag_su2==1
    su2_in=cell(1,3);
    su2_in{1}=zeros(n,n);
    su2_in{2}=zeros(n,n);
    su2_in{3}=zeros(n,n);
    search_val = cell(1,3);
    search_val{1}='';
    search_val{2}='';
    search_val{3}='';
    
    for i=1:n %цикл по заказам
        %     if i ==3
        %         stop=1;
        %     end
        f_srch = id_in_u(i);
        if f_srch == 0
            continue;
        end
        %      if f_srch == 100
        %         stop=1;
        %      end
        id_n = find(id_t==f_srch);
        %cyl=-1;
        search_val{1}='';
        search_val{2}='';
        search_val{3}='';
        k=0;
        
        for l=1:numel(id_n) %цикл по строкам
            ll=id_n(l);
            %if id_in(ll)==f_srch % i
            if strcmp(db_data1{ll,3},'EBC::FX')
                k=k+1;
            end
            if strcmp(db_data1{ll,3},'G03::CU')
                if k==0
                    k=k+1;
                end
                search_val{k} = db_data1{ll,28};
            end
            %end
        end %l
        %если только 2 прохода - 2-й считается как 3-й
        if k==2
            search_val{3}= search_val{2};
            search_val{2} = '';
        end
        for jj=i+1:n %цикл по заказам
            %         if j ==10
            %             stop=1;
            %         end
            f_srch = id_in_u(jj);
            if f_srch == 0
                continue;
            end
            %         if f_srch == 106
            %             stop=1;
            %         end
            id_n = find(id_t==f_srch);
            k=0;
            tmp_l = 0;
            for l=1:numel(id_n) %цикл по строкам
                ll=id_n(l);
                %if id_in(ll)==f_srch % j
                %db_data1{l,3}
                if strcmp(db_data1{ll,3},'EBC::FX')
                    k=k+1;
                end
                if strcmp(db_data1{ll,3},'G03::CU')
                    if k==0
                        k=k+1;
                    end
                    if k==2
                        tmp_l = ll; %сохраним, для еще одной проверки
                    end
                    if ~strcmp(db_data1{ll,28},search_val{k}) && ~strcmp(search_val{k},'')
                        su2_in{k}(i,jj)=10;
                        su2_in{k}(jj,i)=10;
                    end
                end
                %end
            end % l
            %если только 2 прохода - 2-й считается как 3-й
            if k==2
                if tmp_l>0
                    if ~strcmp(db_data1{tmp_l,28},search_val{3}) && ~strcmp(search_val{3},'')
                        su2_in{k}(i,jj)=10;
                        su2_in{k}(jj,i)=10;
                    end
                end
                
                su2_in{3}(i,jj)=su2_in{2}(i,jj);
                su2_in{3}(jj,i)=su2_in{2}(jj,i);
                su2_in{2}(i,jj)=0;
                su2_in{2}(jj,i)=0;
            end
        end % j
    end % i
end
end

