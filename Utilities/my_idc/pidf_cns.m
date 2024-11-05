function  [h,val]=PIDF_cns(K,tol,f_,op_t,conx_);
%  function  [h,val]=PIDF_cns(K,tol,f_,op_t,con_);
%    specifying the PID parameter constraints
%

% K. TSAKALIS, 8/23/96

con_=conx_;h=0*K; val=0;T=f_(1);TF=f_(2);
if TF==0;TF=1;end
if op_t(2)==1;con_(3)=2*tol;end

if con_(4)==0
    if K(1)<=tol,h(1)=-1;val=tol-K(1);          % positivity
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);    % positivity
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);    % positivity
    end
    if K(3)>con_(2)*T*TF & con_(2)>=tol,
        h=[0 0 1]';val=K(3)-con_(2)*T*TF;         % HF bound
    end
    if K(1)-T*K(2)<=tol,
        h=-[1, -T, 0]';val=tol-(K(1)-T*K(2));    % positivity
    end
    if K(3)-T*K(1)+T*T*K(2) <=tol,
        h=-[-T, T*T,1]';
        val=tol-(K(3)-T*K(1)+T*T*K(2));          % positivity
    end
    if K(3)-T*K(2)+T*T*K(1) >= con_(3) & con_(3)>=2*tol,
        h=[T*T,-T,1]';
        val=(K(3)-T*K(2)+T*T*K(1))-con_(3);      % kd bound
    end
else
    if K(1)<=tol,h(1)=-1;val=tol-K(1);
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
    end
    if K(3)>con_(2)*T*TF & con_(2)>=tol,
        h=[0,0,1]';val=K(3)-con_(2)*T*TF;
    end
    if K(3)-T*K(2)+T*T*K(1) >= con_(3) & con_(3)>=2*tol,
        h=[T*T,-T,1]';
        val=(K(3)-T*K(2)+T*T*K(1))-con_(3);
    end
end
