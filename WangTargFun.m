function [SC]=WangTargFun(a1,b1,c1,S)

ns=numel(S);
SC=0;
a_s=sum(a1);
%b_s=sum(b1);
%c_s=sum(c1);
constr=a_s;
e1=zeros(1,ns);
e2=zeros(1,ns);
e3=zeros(1,ns);
phi=zeros(1,ns);
phi2=zeros(1,ns);

for i=1:ns
    j=S(i);
    if i==1
        e1(i)=a1(j);
        e2(i)=e1(i)+b1(j);
        
        if ns==1
            if e2(i)>constr
                e3(i-1)=e2(i);
            elseif i==ns
                e3(i-1)=constr;
            end
            if e3(i-1)>0
                e3(i)=e3(i-1)+c1(j);
                SC=e3(i);
            end
        end
        
    else
        e1(i)=e1(i-1)+a1(j);
        phi(i)=max(e1(i),e2(i-1));
        e2(i)=phi(i)+b1(j);
        constr=constr+c1(S(i-1));
        
        if e3(i-1)==0
            %constr=constr+c1(j);
            if e2(i)>constr
                e3(i-1)=e2(i);
            elseif i==ns
                e3(i-1)=constr;
            end
            if e3(i-1)>0
                e3(i)=e3(i-1)+c1(j);
                SC=e3(i);
            end
        else
            phi2(i)=max(e3(i-1),e2(i));
            e3(i)=phi2(i)+c1(j);
            SC=e3(i);
        end
    end
end
end