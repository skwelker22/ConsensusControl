function [x,y]=digitize(nfig);
%function [x,y]=digitize(nfig);

x=[];y=[];
disp(['Data recording from Fig.  ',num2str(nfig)])
disp('left button: record data point; right button delete last point')
disp('any key to stop data recording')

figure(nfig)
fin=0;

while fin==0
  [xr,yr,b]=ginput(1);
    if b==1
      x=[x;xr];y=[y; yr];
    elseif b==2
      x=x(1:length(x)-1);y=y(1:length(y)-1);
    else
      fin=1;
    end
end
hold on
plot(x,y,'--m')
hold off
