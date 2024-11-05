function  [h,val]=PIDqcns(K,tol,T,c_);
%  function  [h,val]=PIDqcns(K,tol,T,c_);
%    specifying the PID parameter constraints
%

% K. TSAKALIS, 8/11/04

cc=inf;
if c_==1;cc=2*tol;end

h=0*K; val=0;
if K(1)<=tol,h(1)=-1;val=tol-K(1);
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
end


  if K(1)<=tol,h(1)=-1;val=tol-K(1);
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
  end
  if K(3)-T*K(2)+T*T*K(1) >= cc 
    h=[T*T,-T,1]';
    val=(K(3)-T*K(2)+T*T*K(1))-cc;
  end
