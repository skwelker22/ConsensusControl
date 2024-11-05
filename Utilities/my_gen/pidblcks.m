function [ret,x0,str,ts,xts]=pidblcks(t,x,u,flag);
%PIDBLCKS	is the M-file description of the SIMULINK system named PIDBLCKS.
%	The block-diagram can be displayed by typing: PIDBLCKS.
%
%	SYS=PIDBLCKS(T,X,U,FLAG) returns depending on FLAG certain
%	system values given time point, T, current state vector, X,
%	and input vector, U.
%	FLAG is used to indicate the type of output to be returned in SYS.
%
%	Setting FLAG=1 causes PIDBLCKS to return state derivatives, FLAG=2
%	discrete states, FLAG=3 system outputs and FLAG=4 next sample
%	time. For more information and other options see SFUNC.
%
%	Calling PIDBLCKS with a FLAG of zero:
%	[SIZES]=PIDBLCKS([],[],[],0),  returns a vector, SIZES, which
%	contains the sizes of the state vector and other parameters.
%		SIZES(1) number of states
%		SIZES(2) number of discrete states
%		SIZES(3) number of outputs
%		SIZES(4) number of inputs
%		SIZES(5) number of roots (currently unsupported)
%		SIZES(6) direct feedthrough flag
%		SIZES(7) number of sample times
%
%	For the definition of other parameters in SIZES, see SFUNC.
%	See also, TRIM, LINMOD, LINSIM, EULER, RK23, RK45, ADAMS, GEAR.

% Note: This M-file is only used for saving graphical information;
%       after the model is loaded into memory an internal model
%       representation is used.

% the system will take on the name of this mfile:
sys = mfilename;
new_system(sys)
simver(1.3)
if (0 == (nargin + nargout))
     set_param(sys,'Location',[273,85,773,385])
     open_system(sys)
end;
set_param(sys,'algorithm',     'RK-45')
set_param(sys,'Start time',    '0.0')
set_param(sys,'Stop time',     '999999')
set_param(sys,'Min step size', '0.0001')
set_param(sys,'Max step size', '10')
set_param(sys,'Relative error','1e-3')
set_param(sys,'Return vars',   '')

add_block('built-in/Note',[sys,'/','set-point'])
set_param([sys,'/','set-point'],...
		'Font Name','Times New Roman',...
		'Font Angle','italic',...
		'Font Size',8,...
		'position',[135,35,140,40])

add_block('built-in/Note',[sys,'/','output'])
set_param([sys,'/','output'],...
		'Font Angle','italic',...
		'Font Size',8,...
		'position',[135,74,140,79])


%     Subsystem  'standard PID'.

new_system([sys,'/','standard PID'])
set_param([sys,'/','standard PID'],'Location',[4,62,604,488])

add_block('built-in/Inport',[sys,'/','standard PID/System Output'])
set_param([sys,'/','standard PID/System Output'],...
		'Port','2',...
		'position',[30,115,50,135])

add_block('built-in/Gain',[sys,'/','standard PID/Proportional'])
set_param([sys,'/','standard PID/Proportional'],...
		'Gain','kp',...
		'position',[265,329,310,361])

add_block('built-in/Constant',[sys,'/','standard PID/Constant'])
set_param([sys,'/','standard PID/Constant'],...
		'Value','0',...
		'position',[375,301,415,329])

add_block('built-in/Gain',[sys,'/','standard PID/gain'])
set_param([sys,'/','standard PID/gain'],...
		'Gain','DT*kp/ki',...
		'position',[265,159,350,211])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','standard PID/Low-pass filter'])
set_param([sys,'/','standard PID/Low-pass filter'],...
		'Numerator','[DT/(TF+DT),  0]',...
		'Denominator','[1 -TF/(TF+DT)]',...
		'Sample time','DT',...
		'position',[115,34,225,86])

add_block('built-in/Inport',[sys,'/','standard PID/Set point'])
set_param([sys,'/','standard PID/Set point'],...
		'position',[30,40,50,60])

add_block('built-in/Sum',[sys,'/','standard PID/P+I+D'])
set_param([sys,'/','standard PID/P+I+D'],...
		'inputs','++++-',...
		'position',[500,265,520,425])

add_block('built-in/Sum',[sys,'/','standard PID/error'])
set_param([sys,'/','standard PID/error'],...
		'inputs','+-',...
		'position',[80,42,100,78])

add_block('built-in/Outport',[sys,'/','standard PID/Input to system'])
set_param([sys,'/','standard PID/Input to system'],...
		'position',[670,332,685,358])

add_block('built-in/Gain',[sys,'/','standard PID/Proportional3'])
set_param([sys,'/','standard PID/Proportional3'],...
		'Gain','kp*kd/SAM',...
		'position',[130,355,240,395])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','standard PID/Dis. Transfer Fcn1'])
set_param([sys,'/','standard PID/Dis. Transfer Fcn1'],...
		'Numerator','[DT/SAM]',...
		'Denominator','[1 -1+DT/SAM]',...
		'Sample time','DT',...
		'position',[305,411,385,449])


%     Subsystem  'standard PID/integrator1'.

new_system([sys,'/','standard PID/integrator1'])
set_param([sys,'/','standard PID/integrator1'],'Location',[4,42,646,496])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','standard PID/integrator1/delay'])
set_param([sys,'/','standard PID/integrator1/delay'],...
		'orientation',2,...
		'Denominator','[1 0]',...
		'Sample time','DT',...
		'position',[250,39,295,81])

add_block('built-in/Outport',[sys,'/','standard PID/integrator1/out_1'])
set_param([sys,'/','standard PID/integrator1/out_1'],...
		'position',[480,130,500,150])

add_block('built-in/Inport',[sys,'/','standard PID/integrator1/in_1'])
set_param([sys,'/','standard PID/integrator1/in_1'],...
		'position',[15,140,35,160])

add_block('built-in/Sum',[sys,'/','standard PID/integrator1/error1'])
set_param([sys,'/','standard PID/integrator1/error1'],...
		'position',[65,122,85,158])
add_line([sys,'/','standard PID/integrator1'],[90,140;475,140])
add_line([sys,'/','standard PID/integrator1'],[40,150;60,150])
add_line([sys,'/','standard PID/integrator1'],[420,140;420,60;300,60])
add_line([sys,'/','standard PID/integrator1'],[245,60;40,60;40,130;60,130])


%     Finished composite block 'standard PID/integrator1'.

set_param([sys,'/','standard PID/integrator1'],...
		'position',[405,159,435,211])
add_line([sys,'/','standard PID'],[420,315;495,315])
add_line([sys,'/','standard PID'],[355,185;400,185])
add_line([sys,'/','standard PID'],[230,60;250,60;260,185])
add_line([sys,'/','standard PID'],[315,345;495,345])
add_line([sys,'/','standard PID'],[245,375;495,375])
add_line([sys,'/','standard PID'],[390,430;435,430;435,405;495,405])
add_line([sys,'/','standard PID'],[230,60;240,60;240,175;105,175;105,375;125,375])
add_line([sys,'/','standard PID'],[245,375;265,375;265,430;300,430])
add_line([sys,'/','standard PID'],[230,60;240,60;240,345;260,345])
add_line([sys,'/','standard PID'],[55,125;60,125;60,70;75,70])
add_line([sys,'/','standard PID'],[55,50;75,50])
add_line([sys,'/','standard PID'],[105,60;110,60])
add_line([sys,'/','standard PID'],[440,185;460,185;460,285;495,285])
add_line([sys,'/','standard PID'],[525,345;665,345])
set_param([sys,'/','standard PID'],...
		'Mask Display','PID+filter',...
		'Mask Type','Discrete PID ')
set_param([sys,'/','standard PID'],...
		'Mask Dialogue','Discrete PID with lowpass filter|Gain kp - Integral time ki - Derivative time kd|Pseudo-derivative filter time constant|Error lowpass filter Time constant|Sampling Time')
set_param([sys,'/','standard PID'],...
		'Mask Translate','kp=@1(1);ki=@1(2);kd=@1(3);SAM=@2;TF=@3;DT=@4;',...
		'Mask Entries','[.05, 3, 0.5]\/1/60\/0.1\/1/60\/')


%     Finished composite block 'standard PID'.

set_param([sys,'/','standard PID'],...
		'position',[190,34,285,96])


%     Subsystem  'AW-PID '.

new_system([sys,'/','AW-PID '])
set_param([sys,'/','AW-PID '],'Location',[4,62,604,488])

add_block('built-in/Inport',[sys,'/','AW-PID /System Output'])
set_param([sys,'/','AW-PID /System Output'],...
		'Port','2',...
		'position',[30,115,50,135])

add_block('built-in/Gain',[sys,'/','AW-PID /Proportional'])
set_param([sys,'/','AW-PID /Proportional'],...
		'Gain','kp',...
		'position',[265,329,310,361])

add_block('built-in/Constant',[sys,'/','AW-PID /Constant'])
set_param([sys,'/','AW-PID /Constant'],...
		'Value','0',...
		'position',[375,301,415,329])

add_block('built-in/Gain',[sys,'/','AW-PID /gain'])
set_param([sys,'/','AW-PID /gain'],...
		'Gain','DT*kp/ki',...
		'position',[265,159,350,211])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','AW-PID /Low-pass filter'])
set_param([sys,'/','AW-PID /Low-pass filter'],...
		'Numerator','[DT/(TF+DT),  0]',...
		'Denominator','[1,  -TF/(TF+DT)]',...
		'Sample time','DT',...
		'position',[115,34,225,86])

add_block('built-in/Inport',[sys,'/','AW-PID /Set point'])
set_param([sys,'/','AW-PID /Set point'],...
		'position',[30,40,50,60])

add_block('built-in/Sum',[sys,'/','AW-PID /P+I+D'])
set_param([sys,'/','AW-PID /P+I+D'],...
		'inputs','++++-',...
		'position',[500,265,520,425])

add_block('built-in/Sum',[sys,'/','AW-PID /error'])
set_param([sys,'/','AW-PID /error'],...
		'inputs','+-',...
		'position',[80,42,100,78])

add_block('built-in/Outport',[sys,'/','AW-PID /Input to system'])
set_param([sys,'/','AW-PID /Input to system'],...
		'position',[670,332,685,358])

add_block('built-in/Saturation',[sys,'/','AW-PID /Saturation'])
set_param([sys,'/','AW-PID /Saturation'],...
		'Lower Limit','umin',...
		'Upper Limit','umax',...
		'position',[570,335,595,355])

add_block('built-in/Gain',[sys,'/','AW-PID /Proportional3'])
set_param([sys,'/','AW-PID /Proportional3'],...
		'Gain','kp*kd/SAM',...
		'position',[130,355,240,395])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','AW-PID /Dis. Transfer Fcn1'])
set_param([sys,'/','AW-PID /Dis. Transfer Fcn1'],...
		'Numerator','[DT/SAM]',...
		'Denominator','[1 -1+DT/SAM]',...
		'Sample time','DT',...
		'position',[305,411,385,449])


%     Subsystem  'AW-PID /integrator2'.

new_system([sys,'/','AW-PID /integrator2'])
set_param([sys,'/','AW-PID /integrator2'],'Location',[4,42,646,496])

add_block('built-in/Switch',[sys,'/','AW-PID /integrator2/Switch2'])
set_param([sys,'/','AW-PID /integrator2/Switch2'],...
		'position',[390,114,415,166])

add_block('built-in/Switch',[sys,'/','AW-PID /integrator2/Switch'])
set_param([sys,'/','AW-PID /integrator2/Switch'],...
		'position',[220,109,245,161])

add_block('built-in/Sum',[sys,'/','AW-PID /integrator2/error2'])
set_param([sys,'/','AW-PID /integrator2/error2'],...
		'inputs','+-',...
		'position',[155,123,175,147])

add_block('built-in/Sum',[sys,'/','AW-PID /integrator2/error3'])
set_param([sys,'/','AW-PID /integrator2/error3'],...
		'inputs','+-',...
		'position',[315,128,335,152])

add_block('built-in/Discrete Transfer Fcn',[sys,'/','AW-PID /integrator2/delay'])
set_param([sys,'/','AW-PID /integrator2/delay'],...
		'orientation',2,...
		'Denominator','[1 0]',...
		'Sample time','DT',...
		'position',[180,39,225,81])

add_block('built-in/Constant',[sys,'/','AW-PID /integrator2/Constant3'])
set_param([sys,'/','AW-PID /integrator2/Constant3'],...
		'Value','umin',...
		'position',[210,193,285,227])

add_block('built-in/Constant',[sys,'/','AW-PID /integrator2/Constant'])
set_param([sys,'/','AW-PID /integrator2/Constant'],...
		'Value','umax',...
		'position',[45,180,115,220])

add_block('built-in/Outport',[sys,'/','AW-PID /integrator2/out_1'])
set_param([sys,'/','AW-PID /integrator2/out_1'],...
		'position',[480,130,500,150])

add_block('built-in/Inport',[sys,'/','AW-PID /integrator2/in_1'])
set_param([sys,'/','AW-PID /integrator2/in_1'],...
		'position',[15,140,35,160])

add_block('built-in/Sum',[sys,'/','AW-PID /integrator2/error1'])
set_param([sys,'/','AW-PID /integrator2/error1'],...
		'position',[65,122,85,158])
add_line([sys,'/','AW-PID /integrator2'],[40,150;60,150])
add_line([sys,'/','AW-PID /integrator2'],[420,140;475,140])
add_line([sys,'/','AW-PID /integrator2'],[420,140;430,140;430,60;230,60])
add_line([sys,'/','AW-PID /integrator2'],[250,135;280,135;280,100;355,100;355,125;385,125])
add_line([sys,'/','AW-PID /integrator2'],[90,140;100,140;100,120;215,120])
add_line([sys,'/','AW-PID /integrator2'],[120,200;180,200;180,150;215,150])
add_line([sys,'/','AW-PID /integrator2'],[250,135;310,135])
add_line([sys,'/','AW-PID /integrator2'],[290,210;350,210;350,155;385,155])
add_line([sys,'/','AW-PID /integrator2'],[290,210;300,210;300,185;275,185;275,145;310,145])
add_line([sys,'/','AW-PID /integrator2'],[340,140;385,140])
add_line([sys,'/','AW-PID /integrator2'],[90,140;150,140])
add_line([sys,'/','AW-PID /integrator2'],[120,200;130,200;130,130;150,130])
add_line([sys,'/','AW-PID /integrator2'],[180,135;215,135])
add_line([sys,'/','AW-PID /integrator2'],[175,60;40,60;40,130;60,130])
set_param([sys,'/','AW-PID /integrator2'],...
		'Mask Display','DT-I-AW',...
		'Mask Type','Anti windup intgrator',...
		'Mask Dialogue','Integral with anti-windup  |umax|umin|sampling time',...
		'Mask Translate','umax=@1 ; umin =@2 ;DT=@3;')
set_param([sys,'/','AW-PID /integrator2'],...
		'Mask Entries','umax\/umin\/DT\/')


%     Finished composite block 'AW-PID /integrator2'.

set_param([sys,'/','AW-PID /integrator2'],...
		'position',[385,159,440,211])
add_line([sys,'/','AW-PID '],[420,315;495,315])
add_line([sys,'/','AW-PID '],[600,345;665,345])
add_line([sys,'/','AW-PID '],[525,345;565,345])
add_line([sys,'/','AW-PID '],[355,185;380,185])
add_line([sys,'/','AW-PID '],[230,60;250,60;260,185])
add_line([sys,'/','AW-PID '],[315,345;495,345])
add_line([sys,'/','AW-PID '],[245,375;495,375])
add_line([sys,'/','AW-PID '],[390,430;435,430;435,405;495,405])
add_line([sys,'/','AW-PID '],[230,60;240,60;240,175;105,175;105,375;125,375])
add_line([sys,'/','AW-PID '],[245,375;265,375;265,430;300,430])
add_line([sys,'/','AW-PID '],[230,60;240,60;240,345;260,345])
add_line([sys,'/','AW-PID '],[55,125;60,125;60,70;75,70])
add_line([sys,'/','AW-PID '],[55,50;75,50])
add_line([sys,'/','AW-PID '],[105,60;110,60])
add_line([sys,'/','AW-PID '],[445,185;460,185;460,285;495,285])
set_param([sys,'/','AW-PID '],...
		'Mask Display','AW-PID+filter',...
		'Mask Type','Anti-windup PID ')
set_param([sys,'/','AW-PID '],...
		'Mask Dialogue','Discrete PID with Anti-windup and with lowpass filter|Gain kp - Integral time ki - Derivative time kd|Pseudo-derivative filter time constant|Error lowpass filter Time constant|Sampling Time|Saturation Limits [umax umin]')
set_param([sys,'/','AW-PID '],...
		'Mask Translate','kp=@1(1);ki=@1(2);kd=@1(3);SAM=@2;TF=@3;DT=@4;umax=@5(1);umin=@5(2);',...
		'Mask Entries','[.05, 3, 0.5]\/1/60\/0.1\/1/60\/[0.7 -0.3]\/')


%     Finished composite block 'AW-PID '.

set_param([sys,'/','AW-PID '],...
		'position',[335,34,425,96])

drawnow

% Return any arguments.
if (nargin | nargout)
	% Must use feval here to access system in memory
	if (nargin > 3)
		if (flag == 0)
			eval(['[ret,x0,str,ts,xts]=',sys,'(t,x,u,flag);'])
		else
			eval(['ret =', sys,'(t,x,u,flag);'])
		end
	else
		[ret,x0,str,ts,xts] = feval(sys);
	end
else
	drawnow % Flash up the model and execute load callback
end
