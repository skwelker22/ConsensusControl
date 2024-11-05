function w_E=u_d_app(fr,mag,nout);

ord_app=ones(1,nout)*3;
EUNC=mag*ones(1,nout);

    while ord_app(1)>0
       w_E=fitd(log(EUNC'),fr',ord_app,ones(1,nout));
       if flag <=2;pause(2);else;pause;end
%       [AWE,BWE,CWE,DWE]=branch(w_E);
%       [AWEi,BWEi,CWEi,DWEi]=ssinv(AWE,BWE,CWE,DWE);
       ord_appx=input('enter app. order, 0 if done [0]  ');
         if length(ord_appx)<1;ord_app=0;
         else ord_app(1:length(ord_appx))=ord_appx;
         end
       clf
    end
