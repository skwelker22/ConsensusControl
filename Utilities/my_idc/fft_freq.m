function frun=fft_freq(no_pts,step_siz);

%  function frun=fft_freq(no_pts,step_siz);
%  to generate the continuous frequencies in rad/s for fft data 
%  f(1)=f(2)/20, an arbitrary choice

if nargin == 1
    t=no_pts;
    no_pts=length(t); step_siz=mean(diff(t));
end

frun=([1:no_pts]'-1)*2*pi/no_pts/step_siz;frun(1)=frun(2)/20;
