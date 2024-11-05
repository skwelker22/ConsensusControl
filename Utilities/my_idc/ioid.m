% script file ioid: Performs I/O system id for a SISO system.
% Requires the plant input "u" and output "y" to be in the
%    workspace.  
% It is recommended that u&y have radix-2 lengths for
%    faster computations.  To avoid excessive bias,
%    the test input should be selected so that the 
%    system is approximately "at rest" in the beginning 
%    and the end of the id-process.  Also, it is assumed
%    that u&y have been pre-processed to remove constant
%    bias terms (i.e., at rest means (u,y)=(0,0))
% Results: numo,deno: identified system numerator and denominator
%          frw, magno, phaso: frequency response (parametric)
%          fftfr, frpla: system frequency response (fft)
%          fftfr, fftmunc: mult.uncert. bound
%          fftfr, fftmync: feedbk unc. bound

% -------- Initialization

hold off, format short e
stp=input('time step  ');
nn=length(y);frw=logspace(-2,2.5,100)';frg=logspace(-1.5,1.3,60)';
valid=input('identification (1) or validation (0)  ');
nc=input('Model type (0=strictly proper,-1=biproper)  ');
filord=input('filter order  ');filpol=input('filter cutoff  ');
prenum=input('prefilter zeros  ');preden=input('prefilter poles  ');
t=[0:nn-1]'*stp;clg;subplot(121);plot(t,y);subplot(122);plot(t,u);pause

% -------- auxiliary filter definition 

den=poly(-ones(1,filord)*filpol);num=[10];
[f,q,cc,dd]=tf2ss(num,real(den));
cy=eye(length(q),length(q));zq=q*0;n=length(q);

% -------- prefilter definition
nuf=poly([-prenum])*preden(1)/prenum(1); def=poly([-preden]);
   if prenum~=preden
      uf=lsim(nuf,def,u,t); yf=lsim(nuf,def,y,t);
      else, uf=u;yf=y;
   end

% -------- filter states
wu=lsim(f,q,cy,zq,uf,t);wy=lsim(f,q,cy,zq,yf,t);
disp('lsim results')
www=[wu,wy];if nc == -1,www=[www,uf]; end

% -------- parameter estimation
   if valid == 1
     	thx=inv(www'*www)*www'*(yf);thx=real(thx);
      th1=thx(1:n);th2=thx(n+1:2*n);
      th3=-thx(length(thx))*min(0,nc);
   end

% -------- Display Results

disp('Estimated parameters')
[th1',th2',th3]             %rem: lth=-th3/(th2'*inv(f)*q)
[numo,deno]=ss2tf(f'+th2*q',th1+th2*th3,q',th3,1);
[magno,phaso]=bode(numo,deno,frw);
disp('Numerator'),numo
disp('Denominator'),deno
disp('Poles'),z1=roots(deno)
disp('Zeros'),z2=roots(numo)

% -------- Compute estimation error without prefiltering

wu=lsim(f,q,cy,zq,u,t);wy=lsim(f,q,cy,zq,y,t);
err=y-[wu,wy,u]*[th1;th2;th3];

% -------- Frequency domain analysis

disp('fft analysis')
fftu=fft(u);fftu=fftu(1:nn/4);fftu=fftu-u(1)/2-u(nn)/2;
ffty=fft(y);ffty=ffty(1:nn/4);ffty=ffty-y(1)/2-y(nn)/2;
fftfr=[0:nn/4-1]'*2*pi/nn/stp;fftfr(1)=fftfr(2)/10;
frpla=ffty./fftu;    %rem: fft estimate of freq.resp.

disp('uncertainty computations')
ffte=fft(err);ffte=ffte(1:nn/4);ffte=ffte-err(1)/2-err(nn)/2;
frrespo=freqs(numo,deno,fftfr);frrespn=freqs(numo,den,fftfr);
frrespd=freqs(den,deno,fftfr);
fftmunc=abs(ffte./fftu)./abs(frrespn);fftmync=abs(ffte./ffty).*abs(frrespd);
fmunc=smoothin(fftfr,fftmunc,frg);fmync=smoothin(fftfr,fftmync,frg);
% -------- PLOTS

hold off, clg,  plot(t,err), title('prediction error')
pause

clg,subplot(221),loglog(frw,magno),grid,title('frequency response')
subplot(223), semilogx(frw,phaso), grid
subplot(122), z1r=real(z1);z1i=imag(z1);
  if length(z2) > 0
    z2r=real(z2);z2i=imag(z2);
    plot(z1r,z1i,'x',z2r,z2i,'o'),grid,title('pole-zero plots')
  else
    plot(z1r,z1i,'x'), grid,title('pole-zero plots')
  end
pause

clg,loglog(frw,magno,fftfr,abs(frpla))
title('magnitude of freq. resp.'),grid,pause
loglog(frg,(fmunc.^(-1)))
title('estimate of inverse multiplicative unc.(T-bound)'),grid,pause
loglog(frg,abs(fmync.^(-1)))
title('estimate of inverse feedback unc. (S-bound)'),grid,pause




