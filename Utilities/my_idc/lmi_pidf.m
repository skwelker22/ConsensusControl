function [K,Ub]=LMI_PIDF(W,zi,toler,f_,op_t,con_,fr,flag);
% function [K,Ub]=LMI_PIDF(W,zi,toler,f_,op_t,con_,fr,flag);
%  Convex optimization function to perform PID tuning
%  Calls to LMI_upd, PIDf_cns, PIDf_obj
%

% K. TSAKALIS, 8/23/96

if nargin<8;flag =0;end

K=[0;0;0]; Aell=con_(1)*con_(1)*eye(length(K),length(K));
L=[];U=[];err=1;erro=1;kount=0;phi=0;kkount=0;indp=0;
op_tx=op_t;if op_tx(1) < 0;op_tx(1)=1;end

while err>toler
    [h,psi]=pidf_cns(K,1.e-5,f_,op_t,con_);
    if norm(h) == 0 & op_t(1) < 0
        [h,psi]=pidf_obj(W,zi,K,op_t,fr);
    end
    h0=norm(h);
    if h0>0
        % 'constraint-iteration'
        iterty=-1;deepcut=psi;
    else
        % 'objective-iteration'
        iterty=1;
        [h,phi]=pidf_obj(W,zi,K,op_tx,fr);
        erro=sqrt(h'*Aell*h);deepcut=0;
        %     L=[L;phi-erro];U=[U;phi];
        if length(U)<1;U=phi;else;U=min(U,phi);end
    end
    [K,Aell,err_flag]=lmi_upd(K,h,Aell,deepcut,flag);
    kount=kount+1;kkount=kkount+1;
    h0=norm(pidf_cns(K,1.e-5,f_,op_t,con_));
    if op_t(1) < 0
        [hx,h1]=pidf_obj(W,zi,K,op_t,fr);
    else
        h1=0;
    end
    err=(h0)+h1+erro;
    if flag >= 2
        if kkount>=10,kkount=0;
            disp([kount, err, phi,iterty])
        end
    end
    if err_flag<0
        if flag>=0;disp('LMI_PIDF WARNING: Optimization failure');end
        err=0;
    end
end

if err_flag<0
    Ub=inf;
else
    Ub=sqrt(U);
end

