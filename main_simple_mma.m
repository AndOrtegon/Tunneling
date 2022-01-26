clear;
clc;

addpath(genpath('util'))

%% Declaration of global variable
global nelx nely volfrac penal rmin ft E0 Emin KE F U H Hs edofMat freedofs iK jK ml penal_type Hproj Hbeta Heta problem epsilon lambda movp
global method 
    
movp = false ;
method = 'direct' ; 

%% Optimization input
% Problem 'Compliance' or 'Heat'
% problem = 'Compliance' ;
% 
% nelx    = 150;
% nely    = 50;
% volfrac = 0.3;
% penal   = 3;
% rmin    = 2;
% ft      = 2;
% cont    = true;
% penal_type = 'RAMP';
% lambda  = 1.5e-4*nelx*nely ;
% epsilon = 1e-1 ; 

% % Heat problem
problem = 'Heat' ;

nelx    = 800;
nely    = 800;
volfrac = 0.4;
penal   = 3;
rmin    = 2;
ft      = 2;
cont    = true ;
penal_type = 'SIMP';
lambda  = 1.5e-4*nelx*nely ;
epsilon = 2.5e-1 ; 

% Move limit
ml      = 0.2;   %1.0 means essentially no move limits

% Heaviside projection
Hproj   = false;
Hbeta   = 8;
Heta    = 0.5;

%% MATERIAL PROPERTIES
if strcmp(problem,'Compliance')
    E0   = 1;
    Emin = 1e-9;
    nu   = 0.3;
    
    % PREPARE FINITE ELEMENT ANALYSIS
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    
    KE      = 1/(1-nu^2)*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])/24;
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
    iK      = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK      = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
    
    % Define loads and supports: Mitchell Bridge
    F = sparse(2*(nely+1)*(nelx/2+1),1,-1,2*(nely+1)*(nelx+1),1);
    U = zeros(2*(nely+1)*(nelx+1),1);
    fixeddofs = union([2*(nely+1),2*(nely+1)-1],[2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1)-1]);
    
    % Define loads and supports: Half-MMB Beam
    % F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
    % U = zeros(2*(nely+1)*(nelx+1),1);
    % fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
    
    % F = sparse(2*(nely+1),1,-1,2*(nely+1)*(nelx+1),1);
    % U = zeros(2*(nely+1)*(nelx+1),1);
    % fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)-1,2*(nelx+1)*(nely+1)]);
    
    alldofs  = [1:2*(nely+1)*(nelx+1)];
    freedofs = setdiff(alldofs,fixeddofs);
elseif strcmp(problem,'Heat')
    E0 = 1 ;
    Emin = 0.001 ;
    
    KE = [ 2/3 -1/6 -1/3 -1/6
          -1/6  2/3 -1/6 -1/3
          -1/3 -1/6  2/3 -1/6
          -1/6 -1/3 -1/6  2/3];
    
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
    iK      = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
    jK      = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
    
    F         = sparse((nely+1)*(nelx+1),1);
    F(:,1)    = 0.01*1600/(nelx*nely);
    fixeddofs = [nely/2+1-(nely/20):2:nely/2+1+(nely/20)];
    alldofs   = [1:(nely+1)*(nelx+1)];
    freedofs  = setdiff(alldofs,fixeddofs);
end

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k  = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% Initialize Design
rng(1);
% x = rand(nely,nelx);
% x = volfrac*ones(nely,nelx);
x = 0.20*ones(nely,nelx);

TOLF = 1e-5;
already_opt = 0;

%% Optimization phase
% output file
opt_hist = fopen('hist.txt','w');
tunnel_hist = fopen('tunnel_hist.txt','w');
opt_s = 'optimum';
opt_tunel_s = 'optimum_tunnel';

% First optimization
[x_star,f_star] = optimization_phase(x,opt_hist,cont);
after_fix = [num2str(1) '.png'];
fig_name  = [opt_s after_fix];
saveas(gcf,fig_name)
f_old = f_star;

x_star_array=[];
x_star_array = [x_star_array,x_star(:)];

%% Tunneling
for i=2:5
    % Tunneling phase
    [x] = tunneling_phase(x_star_array,f_star,tunnel_hist);
    fig_name = [opt_tunel_s after_fix];
    saveas(gcf,fig_name)
    
    % Optimization phase
    ml = 0.1 ;
    [x_star,f_star] = optimization_phase(x,opt_hist,false);
    if(abs(f_old-f_star)>TOLF*(1+f_star))
        %clear up x_star
        x_star_array = [];
        %append new x_star
        x_star_array = [x_star_array,x_star(:)];
        f_old = f_star;
        after_fix = [num2str(i) '.png'];
        fig_name  = [opt_s after_fix];
        saveas(gcf,fig_name)
        
        disp('tunneling found a lower level!');
        pause(0.5)
    else
        x_star_array = [x_star_array,x_star(:)];
        disp('tunneling found another local minimum, add a pole!');
        pause(0.5)
    end
end

fclose(opt_hist);
fclose(tunnel_hist);
