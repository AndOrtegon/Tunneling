% Load data
clear 
clc

addpath(genpath('util'))

global ml movp method cntp tunneling lambda lambdaM0 s_rand plotting
% load('MBridge_300x100.mat') ;
% load('Initial Conditions/MBridge_150x50.mat') ;
% load('2022-01-12/Optimum 150x50.mat') ;
% method = 'direct' ;
load('Initial Conditions/Heat_300x300_cft.mat') ;
method = 'direct' ;
tunneling = 'gaussian' ;
movp = false ;
cntp = false ;
lambda = 1e1 ;
lambdaM0 = 0e-2 ;
s_rand = false ;
plotting = false ;
rng(2)  

%% Tunneling
opt_hist = fopen('hist.txt','w');
tunnel_hist = fopen('tunnel_hist.txt','w');

% xs = reshape(x_star_v(:,:,1),nely*nelx,[]); % Good minimizers array
% fs = c(1) ; % Good min object values
% xp = reshape(x_star_v(:,:,1:9),nely*nelx,9); % Bad minimizers array; % Bad minimizers array
xs = x_star(:); % Good minimizers array
fs = f_star ;
xp = [] ;

for i=2:2
    % Tunneling phase
    [x, xp, headers, hist] = tunneling_phase(xs,fs,xp,tunnel_hist);
    eval(['hist_' int2str(i) ' = hist;']) ;
    fig_name = [opt_tunel_s after_fix];
    saveas(gcf,fig_name)
    
    % Optimization phase
    ml = 0.1 ;
    [x_star,f_star] = optimization_phase(x,opt_hist,false);
    
    if abs(f_old-f_star) < 1e2*TOLF
        % Found a minimum with the same level
        % Update xs array
        xs = [xs, x_star(:)] ;
        fs = [fs, f_star] ;
        
        if f_star < min(fs) ; f_old = f_star ; end
            
        after_fix = [num2str(i) '.png'];
        fig_name  = [opt_s after_fix];
        saveas(gcf,fig_name)
        
        disp('Tunneling found an equivalent minimum');
        pause(0.1)
    elseif 1e2*TOLF < (f_old-f_star)
        % Found a better minimum
        xp = [xp,xs] ;
        xs = x_star(:);
        fs = f_star ;
        
        after_fix = [num2str(i) '.png'];
        fig_name  = [opt_s after_fix];
        saveas(gcf,fig_name)
        
        disp('Tunneling found better local minimum!');
        pause(0.1)
    else 
        xp = [xp,xs] ;
        xs = x_star(:);
        fs = f_star ;
        
        disp('This tunneling search failed!');
        pause(0.1)
    end
end

fclose(opt_hist);
fclose(tunnel_hist);
