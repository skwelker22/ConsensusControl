%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%add libraries for controller design
addpath(genpath('/Users/swelker/Desktop/swelker/ConsensusControl'));

%% driver for controller
plotFlag=false;

%2d case
nDims=2;
[Tquad_ss2d, Kquad2d]=CreateQuadController(nDims,plotFlag);

%3d case
nDims=3;
[Tquad_ss3d, Kquad3d]=CreateQuadController(nDims,plotFlag);

% %error case
% nDims=4;
% [~,~]=CreateQuadController(nDims);

