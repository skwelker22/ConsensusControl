function [A,B,C,D]=flt_bank(N,min_f,max_f,flag);
%function [A,B,C,D]=flt_bank(N,min_f,max_f,flag);
% creates a N+1-filter bank decomposition of 1, 
% with corner frequencies log-distributed between min_f and max_f

if nargin<4;flag=0;end

f_1=log10(min_f);f_2=log10(max_f);
fn=logspace(f_1,f_2,N);
[A,B,C,D]=tf2ss(fn(1),[1 fn(1)]);

if N>1
   for i=1:N-1
      [a,b,c,d]=tf2ss([fn(i+1)-fn(i) 0],[1 fn(i+1)+fn(i) fn(i+1)*fn(i)]);
      [A,B,C,D]=append(A,B,C,D,a,b,c,d);
   end
end
[a,b,c,d]=tf2ss([1 0],[1 fn(max(N,1))]);
[A,B,C,D]=append(A,B,C,D,a,b,c,d);
B=B*ones(max(N,1)+1,1);D=D*ones(max(N,1)+1,1);
if N==0
   C=[1 1]*C;D=[1 1]*D;
end

A_=mksys(A,B,C,D,'ss');
[A,B,C,D]=h_sysred(A_,[],[],flag,0,[3,1e-5]);

if flag>1
   fr=logspace(f_1-1,f_2+1,100)';
   m=[bode(A,B,C,D,1,fr)];
   loglog(fr,m)
end
