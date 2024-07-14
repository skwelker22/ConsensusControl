%% START
clear all; close all; home;

%format
format long g;
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%%
A = [1, 2; 0 1]
N = 3;
IN = eye(N);
AHat = kron(IN, A)

%% Graph example from Lewis: coop ctrl-opt and adaptive
nStates = 2; %n
nRelativeMeas = 3; %m
A = [0, 0, 1, 0, 0, 0;
     1, 0, 0, 0, 0, 1;
     1, 1, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 1, 0];

D = diag(sum(A'));
L_Lewis = (D - A);
L_Fax = inv(D)*(D-A); %normalization to not add gain

G = eye(size(L_Fax)) - L_Fax; % weighted adjacency matrix

lambdaL = eig(L_Fax);

[G_V, G_lambda, G_W] = eig(G);
er = G_V(:,1); el = G_W(:,1);
er_normalize = er/norm(er'*el);

E = er_normalize * el';

%lemma 6.1
G2 = E + (G-E)^2;

Ln = kron(L_Fax, eye(nRelativeMeas));
Lhat = kron(eye(nRelativeMeas), L_Fax);

%shur transformation
[T, U] = schur(L_Fax);

if(0)
    %plots
    %all eigen values of L lie in a disk of radius 1 centered at 1 + j0
    %this is denoted the "Perron Disk"
    tt = [0:0.01:1];
    figure('Name','Eigen Values of L');
    plot( real(lambdaL), imag(lambdaL), 'bx' ); hold on;
    plot( sin(2*pi.*tt)+1, cos(2*pi.*tt), 'k' );
end