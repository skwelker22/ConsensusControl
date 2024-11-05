function cl_spec=q_spec(frun,MUNC,S_roll,RISK,AGRS,T_roll,figno);
% function cl_spec=q_spec(frun,MUNC,S_roll,RISK,AGRS,T_roll);
% quick spec
if nargin<3;S_roll=2;end
if nargin<4;RISK=1;end
if nargin<5;AGRS=.9;end
if nargin<6;T_roll=3;end
if nargin<7;figno=6;end

cl_spec=bw_app(frun,MUNC,S_roll,RISK,AGRS,T_roll);
S_roll=cl_spec(4);RISK=cl_spec(5);AGRS=cl_spec(6);T_roll=cl_spec(2);

    sp_done=0;
      while sp_done==0
        disp('q_spec: RISK,AGRS,S_roll, T_roll defaults:')
        disp([RISK,AGRS,S_roll,T_roll])
        temp=input('q_spec in: Give RISK,AGRS,S_roll,T_roll ?  ');
          if length(temp)==1; 
               if temp ~=0;RISK=temp;end
            elseif length(temp)==2;
               if temp(1)~=0;RISK=temp(1);end;if temp(2)~=0;AGRS=temp(2);end;
            elseif length(temp)==3;
               if temp(1)~=0;RISK=temp(1);end;if temp(2)~=0;AGRS=temp(2);end;
               if temp(3)~=0;S_roll=temp(3);end
            elseif length(temp)==4;
               if temp(1)~=0;RISK=temp(1);end;if temp(2)~=0;AGRS=temp(2);end;
               if temp(3)~=0;S_roll=temp(3);end;if temp(4)~=0;T_roll=temp(4);end;
          end
        cl_spec=bw_app(frun,MUNC,S_roll,RISK,AGRS,-T_roll);
        disp('q_spec: cl_spec result')
        disp(cl_spec)
            Tsens=awgt_sel(frun,cl_spec(1),cl_spec(2),0);
            Sens=awgt_sel(frun,cl_spec(3),cl_spec(4),-1);
            Tsens=Tsens';Sens=Sens';
    figure(figno);clf
    loglog(frun,Tsens,frun,Sens,frun,[1.0./MUNC]);grid
    title('Weights and preliminary uncertainties')

        sp_done=input('q_spec: Done ?  (0=no, [1])  ');
        if isempty(sp_done);sp_done=1;end
      end
    disp('q_spec: Suggested [T-BW, T-ROLLOFF, S-BW, S-ROLLOFF, RISK, AGRS]')
    disp(cl_spec)
    temp_s=input('q_spec: Enter new specs if desired ?  ');
      if length(temp_s)>0;
        cl_spec(1:length(temp_s))=temp_s;
      end
      if cl_spec(3)==0
        temp_s=sbw_rule(cl_spec);
        cl_spec=[cl_spec(1:2),temp_s,cl_spec(4:length(cl_spec))];
      end

R_A=inv_spec(cl_spec,frun,MUNC);
if length(cl_spec)<6
   cl_spec(5:6)=R_A;
else
   if abs(log(R_A(1)/cl_spec(5)))+abs(R_A(2)-cl_spec(6))>0.002
     cl_spec(5:6)=R_A;
   end
end
