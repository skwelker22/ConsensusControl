function [zsm,zstd]=smoothin(f,z,w,squares);
% [zsm,zstd]=smoothin(f,z,w);
% f: freq in, z: fr.resp. in, w: desired freqs
% zsm: averaged z at w, zstd: st. deviation at w

if nargin<4;squares=0;end
zsm=0*(w);zstd=0*(w);
if squares ~= 0; z=z.*z; end
if f(1) > 0; f=[0;f]; z=[z(1);z];end
f=[f;10*max(max(f),max(w))]; z=[z;0];
WL=w(1:length(w)-1);WH=w(2:length(w));
WM=[0;(WH+WL)/2;max(f)];
for i=1:length(w)
    wlow=WM(i);whi=WM(i+1);
    ind=find(f <= whi & f >= wlow);
    if length(ind) < 2
        ind=[max(find(f < w(i)));min(find(f > w(i)))];
        F=f(ind);Z=z(ind);
        zsm(i) = interp1(F,Z,w(i));
    else
        deltaf=max(f(ind))-min(f(ind));
        df=diff(f(ind)); ztem=z(ind); n=length(ind);
        sz=sum((ztem(1:n-1)+ztem(2:n)).*df/2);
        zsm(i)=sz/deltaf;
        
    end

end

if squares ~=0; zsm=sqrt(zsm); end
