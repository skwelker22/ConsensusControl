function [zsm,zstd]=smoothin(f,z,w);
% [zsm,zstd]=smoothin(f,z,w);
% f: freq in, z: fr.resp. in, w: desired freqs
% zsm: averaged z at w, zstd: st. deviation at w

zsm=0*(w);zstd=0*(w);
zind=find(abs(z)<1.e-8);
z(zind)=0*z(zind)+1.e-8;tem=log(z);
for i=1:length(w)
  k=max(i-1,1);l=min(i+1,length(w));
  ind=find(f <= w(l) & f >= w(k));
    if length(ind) == 0
      ind=[max(find(f <= w(i)));min(find(f >= w(i)))];
    end
  zsm(i)=mean(tem(ind));zstd(i)=std(tem(ind));
end
zsm=exp(zsm);zstd=exp(zstd);
