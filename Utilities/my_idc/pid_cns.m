function  [h,val]=PID_cns(K,tol,T,pos,B,BDZ);
%  function  [h,val]=PID_cns(K,tol,T,pos,B,BDZ);
%    specifying the PID parameter constraints
%

h=0*K; val=0;

if pos==0
  if K(1)<=tol,h(1)=-1;val=tol-K(1);
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
    elseif K(3)>B,h(3)=1;val=K(3)-B;
%    elseif K(2)/K(1) < BDZ;  if BDZ>=0 h(2)=-1;h(1)=BDZ;
%           val=abs(BDZ*K(1)-K(2))/sqrt(1+BDZ*BDZ);end
  end
  if K(1)-T*K(2)<=tol,
       h=-[1, -T, 0]';
       val=tol-(K(1)-T*K(2));
  end
  if K(3)-T*K(2)+T*T*K(1) <=tol, 
       h=-[T*T,-T,1]';
       val=tol-(K(3)-T*K(2)+T*T*K(1));
  end



else
  if K(1)<=tol,h(1)=-1;val=tol-K(1);
    elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
    elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
    elseif K(3)>B,h(3)=1;val=K(3)-B;
  end
end

