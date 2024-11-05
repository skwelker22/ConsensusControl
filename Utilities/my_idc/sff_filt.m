function [A,B,C,D]=sff_filt(ny,wp,wz);
% [A,B,C,D]=sff_filt(ny,wp,wz);
% simple feedforward filter


[A,B,C,D]=tf2ss([1/wz(1),1],[1/wp(1),1]);
np=length(wp);nz=length(wz);
for i=2:ny
    den=[1/wp(min(i,np)),1];num=[1/wz(min(i,nz)),1];
    if num(1) == den(1)
        A1=[];B1=[];C1=[];D1=1;
    else
        [A1,B1,C1,D1]=tf2ss(num,den);
    end
    
    [A,B,C,D]=append(A,B,C,D,A1,B1,C1,D1);

end
