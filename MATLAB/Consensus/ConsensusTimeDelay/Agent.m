%This class is to house the dynamics for the double integrator agent,
%Made this a generic class that can be modified for different types of
%dynamics
classdef Agent < handle

    %agent class public properties
    properties
        states=[];
        stateHistory=[];
        controlHistory=[];
        control=[];
        nDims;
        position=0.0;
        velocity=0.0;
        dT=0.0;
    end
    
    %dynamics should be private properties of the class,
    %unable to change from the object interface
    properties(Access=private)
        % Ac=[0, 1; 0, 0]; 
        Ac=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
        Bc=[0, 0; 0, 0; 1, 0; 0, 1];
        Ad=[];
        Bd=[];
    end

    methods
        %class contstructor
        function obj=Agent(x0,nDims,dt)
            %initialize states and time
            obj.nDims=nDims; obj.dT=dt;
            obj.states=x0; obj.stateHistory=x0;
            obj.position=x0(1:nDims); obj.velocity=x0(nDims+1:end);

            %discretize state space w.r.t. dt, use zoh
            [obj.Ad, obj.Bd]=c2d(obj.Ac,obj.Bc,dt);
        end
        
        %set method for control input
        function []=set.control(obj,value)
            obj.control=value;
        end
        
        %function to propogate teh dynamics forward in time by dT
        function obj=PropagateDynamics(obj)
            %propogate the dynamics according to double integrator
            obj.states=obj.Ad * obj.states + obj.Bd * obj.control;
            %update position and velocity members for readability
            obj.position=obj.states(1:obj.nDims); obj.velocity=obj.states(obj.nDims+1:end);
            %keep track of history for analysis
            obj.stateHistory=[obj.stateHistory,obj.states];
            obj.controlHistory=[obj.controlHistory;obj.control];
        end
        
        %method to get position at a given delay sample
        function delay_position=DelayPosition(obj,tdel)
            delay_position=obj.stateHistory(1:obj.nDims,tdel);
        end

    end

end