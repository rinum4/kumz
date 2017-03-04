function [LocDom3] = KUMZ_get_domP3su(i,j,S_full,Un, p,n,m,ma,d,su,end_time_ex, r)
% получим функцию доминирования
%важно! в данной функции сменилась модель хранения su

%LocDom2 = [];
LocDom3 =[];

%for i=1:numel(Un)
%     parfor j=1:numel(Un) %i+1  %вариант 3: по уникальным но по всем индексам
%         %[U(i)  U(j)]
%         if JeongCheckPropAll2_o4(S_full,Un,Un(i),Un(j), p,n,m,ma,d,su,end_time_ex)==1
%             %LocDom2(k2,:) = [U(j) U(i)];
%             %k2=k2+1;
%             tmp = [Un(j) Un(i)];
%             LocDom2 = [LocDom2; tmp];
%         end
        for l=1:numel(Un) %j+1  %parfor бесмысленнен
            %[Un(i)  Un(j) ,Un(l)]
            if JeongCheckPropAll3_o5(S_full,Un,Un(i),Un(j),Un(l), p,n,m,ma,d,su, end_time_ex,0, r)==1
                %LocDom3(k3,:) = [U(j) U(i) U(l)];
                %k3=k3+1;
                tmp = [Un(j) Un(i) Un(l)];
                LocDom3 = [LocDom3; tmp];
            end
        end
%     end
%end

end

