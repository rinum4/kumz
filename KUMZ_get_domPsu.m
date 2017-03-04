function [LocDom2, LocDom3] = KUMZ_get_domPsu(S_full,Un, p,n,m,ma,d,su,end_time_ex, r)
% получим функцию доминирования
%важно! в данной функции сменилась модель хранения su

LocDom2 = [];
LocDom3 =[];

parfor i=1:numel(Un) %parfor
    [tmp1, tmp2]=KUMZ_get_domP2su(i,S_full,Un, p,n,m,ma,d,su,end_time_ex, r);
    LocDom2 = [LocDom2;tmp1];
    LocDom3 = [LocDom3;tmp2];
end

end

