function [Ar,Br,Cr,Dr]=sys_red(A,B,C,D,def_flag,red_ind,frw,bslev);
%function [Ar,Br,Cr,Dr]=sys_red(A,B,C,D,def_flag,red_ind,frw,bslev);
% General reduction by Schur's method.
%    Slow-fast decompositions are used to remove too fast
%    and too slow dynamics. frw is a frequency vector for
%    computation and display of sensitivity singular values.
%    flags: def_flag~= 0, use defaults
%           red_ind: [n1,n2,n3], n1=hi-freq; n2=low-freq
%            (A should be nonsingular); n3=1 (schmr), 2 (bstschmr)
%            n1,n2<0 uses defaults 400,1/10
%            bslev=[bst/sch.cutoff, shift.floor, bst biprop.threshold]
%                   shift.floor<0 => no shifting

if nargin < 5;def_flag=1;end
if nargin < 6; red_ind=[0 0 2];end
if nargin < 7;frw=0;end
if nargin <8; bslev=[-50,1.e-2,1.e-2];end
if length(bslev)<2;bslev=[bslev,1.e-2];end
if length(bslev)<3;bslev=[bslev,1.e-2];end


if frw ~=0 & def_flag>2
magu=sv3_5(A,B,C,D,1,frw);
end

EAC=eig(A);
TZAC=tzero(A,B,C,D);
  if def_flag>4
    disp('Maximum real parts of poles and zeros')
    disp([max(real(EAC)),max(real(TZAC))])
  end
    Ar=A;Br=B;Cr=C;Dr=D;

%  -------SLOW-FAST DECOMPOSITION (HI-FREQ)------------
if red_ind(1) ~=0
  EAC=eig(Ar);
  if def_flag > 1
    CUTfr=input('cutoff freq  [400]  ');
  else
    if red_ind(1)>0,CUTfr=red_ind(1);else,CUTfr=[];end
  end

  if length(CUTfr)<1;CUTfr=400;end
  CUT=length(find(abs(EAC)<=CUTfr))
  if CUT == 0; CUT=length(Ar);end
  if CUT < length(Ar)
    [Ar,Br,Cr,Dr,AH,BH,CH,DH] = slowfast(Ar,Br,Cr,Dr,CUT);
    disp('Including hf-DC-component in D')
    DH-CH*inv(AH)*BH;
        Dr=Dr+DH-CH*inv(AH)*BH;
  end
  if frw ~=0 & def_flag >2
    magur=sv3_5(Ar,Br,Cr,Dr,1,frw);
    loglog(frw,magu,frw,magur,'--')
    title('System singular values after fast modes reduction'),pause
  end
end

%  -------SLOW-FAST DECOMPOSITION (LOW-FREQ)------------
if red_ind(2) ~= 0
  disp('Low freq. reduction')
  EAC=abs(eig(Ar));
  min(EAC)
  if def_flag >1
    itemp=input('Reject low frequencies (0=no)? [1] ');
    if isempty(itemp);itemp=1;end
  else
    itemp=1;
  end

  if itemp ~= 0
    if def_flag >2
      disp('Magnitude of System  Poles')
      sort(EAC)'
    end
      if def_flag >1
        CUTfr=input('cutoff freq [.1]  ');
      else
        if red_ind(2)>0;CUTfr=red_ind(2);else;CUTfr=[];end
      end
    if length(CUTfr)<1;CUTfr=.1;end
    CUT=length(find(EAC<=CUTfr))
    if CUT==0;CUT=length(Ar);end
    if CUT < length(Ar)
      [AH,BH,CH,DH,Ar,Br,Cr,Dr] = slowfast(Ar,Br,Cr,Dr,CUT);
      disp('Including hf-component in D')
      Dr=Dr+DH;
    end
    if frw~=0 & def_flag>2
      magur=sv3_5(Ar,Br,Cr,Dr,1,frw);
      loglog(frw,magu,frw,magur,'--')
      title('System singular values after slow modes reduction'),pause
    end
  end
end

%  -----SCHUR MODEL REDUCTION---------------
if red_ind(3) >0
  if def_flag >2
    itemp=input('System reduction (0=no)? [1] ');
    if isempty(itemp);itemp=1;end
  else
    itemp=1;
  end

  if itemp ~=0
    xdone=0;
    while xdone==0;
      Ar1=Ar;
      if def_flag >1
        met=input('mod-ord red method: 1=schmr, 2=bstschmr ');
      else
        met=red_ind(3);
      end

    emax=max(real(eig(Ar)));
      if def_flag >1
        disp(['max real system eval =  ',num2str(emax)]);
        delred=input('give delta-shift for mod-ord red [emax+.01] ');
      else
        delred=[];
      end
      if length(delred)<1 ; 
        delred=max(emax+bslev(2),0);
        if bslev(2)<0
           delred=0;
        end
      end
      if delred ~=0
        Ar1=Ar-delred*eye(length(Ar),length(Ar));
      end
      if met == 2
        if def_flag >1
          thrbst=input('bstschmr threshold   ');
        else
          thrbst=[];
        end
        if length(thrbst)<1,thrbst=bslev(3);end
        [UD,SD,VD]=svd(Dr);
        SD=diag(max(diag(SD),thrbst));Drd=UD*SD*VD';
        if def_flag >1
         [Arx,Brx,Crx,Drx,tot,hsv]=bstschmr(Ar1,Br,Cr,Drd,3);
        else
         [Arx,Brx,Crx,Drx,tot,hsv]=bstschmr(Ar1,Br,Cr,Drd,2,bslev(1));
        end
      else
          if def_flag >1 | bslev(1)<0
           [Arx,Brx,Crx,Drx,tot,hsv]=schmr(Ar1,Br,Cr,Dr,3);
          else
           [Arx,Brx,Crx,Drx,tot,hsv]=schmr(Ar1,Br,Cr,Dr,2,bslev(1));
          end
      end
      Arx=Arx+delred*eye(length(Arx),length(Arx));
      if frw~=0 & def_flag>2
        magur=sv3_5(Arx,Brx,Crx,Drx,1,frw);
        if length(magur)==length(magu) 
         loglog(frw,magu,frw,magur,'--')
         title('System singular values after reduction'),pause
        end
      end
      if def_flag >1
        xdone=input('Done with Schur reduction? (0=no) [1]  ');
        if isempty(xdone);xdone=1;end
      else
        xdone=1;
      end
    end
    Ar=Arx;Br=Brx;Cr=Crx;
    if met ~=2, Dr=Drx;end 
  end
end
if norm(Br-0*Br)<1.e-9 | norm(Cr-0*Cr)<1.e-9,Ar=[];Br=[];Cr=[];end
if nargout == 1, Ar=ss(Ar,Br,Cr,Dr); end

