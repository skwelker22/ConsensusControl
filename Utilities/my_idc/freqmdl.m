function [fftfr,Gest]=FREQmdl(y,u,stp,M,nover,fpts);
%function [fftfr,Gest]=FREQmdl(y,u,stp,[M,nover,fpts]);
%  Frequency domain id of a MISO response specified in [u,y], via FFTs
%    u,y=I/O pair, stp=sampling time
%    M=window fraction (total points/M); default=8
%    nover=overlap fraction (window/nover; for M>=4); default=4
%    fpts=high frequency limit (window/fpts); default=4
%

if nargin<4;M=0;end
if nargin<5;nover=0;end
if nargin<6;fpts=0;end
if M == 0;M=8;end
if nover==0;nover=4;end
if fpts==0;fpts=4;end

[nn,ninp]=size(u);[nn1,noutp]=size(y);
     win_pts=fix(nn/M);over_pts=fix(win_pts/nover);
     fr_pts=fix(win_pts/fpts);

if noutp ~= 1; disp('frdommodel requires one output');return;end
if nn ~= nn1;
disp('freqmdl error: lengths of y and u must be equal ');
return;end
fftfr=([1:fr_pts]'-.5)*2*pi/win_pts/stp;
wind=hanning(win_pts);wind2=hanning(fix(nn/2));
tempu=ones(ninp,1);
Cyu=zeros(fr_pts,ninp);Cuu=zeros(ninp*fr_pts,ninp);
nor_fac=sqrt(trace(u'*u));

if M >=1
   fy=fft(y)/nor_fac;fu=fft(u)/nor_fac;
   for i=1:fr_pts
LL=conj(fu((i-1)*M+1:i*M,:));
Cyu(i,:)=Cyu(i,:)+sum((fy((i-1)*M+1:i*M)*tempu').*LL)*10;
Cuu(ninp*(i-1)+1:ninp*i,:)=Cuu(ninp*(i-1)+1:ninp*i,:)+...
    LL'*LL*10;
   end
end

disp('1st pass')

if M>=2
for i=0:4;
   tem_range=i*fix(nn/8)+1:i*fix(nn/8)+fix(nn/2);
   uhan=u(tem_range,:).*(wind2*tempu');
   yhan=y(tem_range,:).*(wind2);
   fy=fft(yhan)/nor_fac;fu=fft(uhan)/nor_fac;
     for j=1:fr_pts
       LL=conj(fu((j-1)*M/2+1:j*M/2,:));
       Cyu(j,:)=Cyu(j,:)+sum((fy((j-1)*M/2+1:j*M/2)*tempu').*LL);
       Cuu(ninp*(j-1)+1:ninp*j,:)=Cuu(ninp*(j-1)+1:ninp*j,:)+...
             LL'*LL;
     end
end
end
disp('2nd pass')


if M>=4

for i=0:nover*(M-1);
   tem_range=i*over_pts+1:i*over_pts+win_pts;
   uhan=u(tem_range,:).*(wind*tempu');
   yhan=y(tem_range,:).*(wind);
   fy=fft(yhan);fu=fft(uhan);
   fy=fy(1:fr_pts)/nor_fac;fu=fu(1:fr_pts,:)/nor_fac;
%  size(fy),size(fu),size(Cyu),size(tempu)
   Cyu=Cyu+(fy*tempu').*conj(fu);
     for j=1:fr_pts
      LL=conj(fu(j,:));
      Cuu(ninp*(j-1)+1:ninp*j,:)=Cuu(ninp*(j-1)+1:ninp*j,:)+...
          LL'*LL;
     end
end
end

disp('3rd pass')

Gest=zeros(fr_pts,ninp);
for i=1:fr_pts
Gest(i,:)=Cyu(i,:)*inv(Cuu((i-1)*ninp+1:i*ninp,:));
end

