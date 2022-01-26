function [x, xp, headers, hist] = tunneling_phase(xs,fs,xp,fileID)
global nelx nely volfrac ml E0 epsilon lambda movp cntp tunneling lambdaM0 s_rand

% It also seems the Heaviside projection was not done in the tunneling
% phase. (why?)

ndv = nelx*nely;
ncons = 1;
ml = 0.2; % move limit

% penal_local = penal;
% Tunneling parameters
beta    = 0.5;
lambdaM = 0 ;
cref    = 0 ;
xm      = ones(nely,nelx) ;

maxoutit = 300;
max_ts   = 5;
counter  = 1 ;
hist     = zeros(max_ts*maxoutit,10) ;
headers  = ["Tunneling step" , "Loop" , "Objective" , ...
    "Constraint" , "KKT-norm" , "Relative change" , ...
    "Tunneling function value" , "Step" , "Distance to previous solution" , ...
    "Movable pole state" ] ;

for tunneling_step = 1:max_ts
    disp('Starting tunneling phase') ;

    x_star = reshape(xs(:,end),nely,nelx);
    
    if s_rand 
        x = randomDirection(x_star,epsilon) ;
    else
        x = volfrac*ones(nely,nelx) ;
    end
    
    xPhys = filterDensities(x) ;
    % xProj = projectDensites(xPhys);
    xProj = xPhys;
    
    % Poles
    pole = [xp,xs] ;
    poleS = lambda*ones(1,size(pole,2)) ;
    
    % INITIALIZE MMA PARAMETERS    
    a0 = 1;
    ai = zeros(ncons,1);
    ci = 1000*ones(ncons,1);
    d  = 0;
    
    xmin = zeros(ndv,1);
    xmax = ones(ndv,1);
    
    low = ones(ndv,1);
    upp = ones(ndv,1);
    
    xold1=x(:);
    xold2=x(:);
    
    kkttol   = 1e-4;
    obj_tol  = 1e-4;
    kktnorm  = kkttol+10;
    old_obj = 0;
    rel_obj_change = 10*obj_tol;
    
    loop = 0;
    while loop < maxoutit
        % Optimality criteria
        if ( rel_obj_change < obj_tol || kktnorm < kkttol ) && (c-fs) > 0
            if movp
                % Convergió
                if lambdaM == 0
                    % Si no había polo movil, encender
                    if 1e-6 < norm(x(:)-xold1(:))
%                         lambdaM = 0.25*dot(x(:)-xold1(:),x(:)-xold1(:)) ;
                        lambdaM = lambdaM0 ;
%                         xm      = xold1(:) ;
                        xm      = x(:) ;
                        x = randomDirection(x,0.1) ;
                    elseif 1e-6 < norm(x(:)-xold2(:))
%                         lambdaM = 0.25*dot(x(:)-xold2(:),x(:)-xold2(:)) ;
                        lambdaM = lambdaM0 ;
                        xm      = x(:) ;
                        x = randomDirection(x,0.1) ;
%                         xm      = xold2(:) ;
                    else
                        break ;
                    end
                    cref    = c ;
                else
                    % Si ya teniamos polo movil, aumentar la fuerza
                    lambdaM = 1.25*lambdaM ;
                end
            else
                % Stopping criteria if movable pole is not used in simulation
                print('Tolerance satisfied') ;
                break;
            end
        elseif movp && 0 < loop && lambdaM ~= 0
            % Else, if using mov pole and pole activated
            if c < cref %|| sign(dot(x(:)-xold1(:),df0dx)) == sign(dot(x(:)-xold1(:),df0dxNP)) ...
                 % || ... 0 < dot(df0dxNP,(x(:)-xm(:)))
                % Conditions to deactivate the movable pole
                lambdaM = 0 ;
            elseif 0.005*nelx*nely < dot(x(:)-xm(:),x(:)-xm(:))
                % Move pole if too far
                xm = 0.75*x(:) - 0.25*xm ;
            end
        end
        
        %% Analize system
        [c, vol, dc, dv] = analyze(xPhys,xProj) ;
        
        xval = x(:);
        xmin = max(0,xval-ml);
        xmax = min(1,xval+ml);
        
        fval  = vol/(volfrac*nely*nelx) - 1;
        dfdx  = dv(:)'/(volfrac*nely*nelx);
        
        %% Stoping criteria with movable pole
        if (rel_obj_change < obj_tol || kktnorm < kkttol) && ~movp && (c-fs) < 0
            % Si es menor 0, termine
            break ;
        end
        
        %% Compute poles and penalizing terms
        % Use a central pole instead of a fixed one
        if cntp
            x_c = beta*x_star(:) + (1-beta)*x(:);
            lambda = 0;
            pole(:,end) = x_c ;
            poleS(end) = lambda;
        end
        
        % Compute product and sum terms of poles
        prod_fix_p = 1;
        sum_fix_p = 0;
        
        for iP = 1:size(pole,2)
            difP = x(:)-pole(:,iP) ;
            dotP = dot(difP,difP) ;
            if strcmp(tunneling,'classical')     
                prod_fix_p = prod_fix_p * dotP^poleS(iP);
                sum_fix_p  = sum_fix_p  + poleS(iP)*difP/dotP^(poleS(iP)+1);
            elseif strcmp(tunneling,'exponential') 
                prod_fix_p = prod_fix_p * exp(poleS(iP)./dotP);
                sum_fix_p  = sum_fix_p  + poleS(iP)*difP/dotP;
            elseif strcmp(tunneling,'gaussian') 
                prod_fix_p = prod_fix_p * ( exp(-dotP/poleS(iP)) + 1 ) ;
                sum_fix_p  = sum_fix_p  + ( difP*exp(-dotP/poleS(iP))/poleS(iP) ) / ...
                    ( exp(-poleS(iP)*dotP) + 1 ) ;
            end
        end
        
        % Movable pole
        prod_mov_p = 1 ;
        sum_mov_p  = 0 ;
        
        if movp
            difM = x(:)-xm(:) ;
            dotM = dot(difM,difM) ;
            if strcmp(tunneling,'classical') 
                prod_mov_p = 1/dotM^lambdaM;
                sum_mov_p  = lambdaM*difM/dotM^(lambdaM+1); 
            elseif strcmp(tunneling,'exponential')         
                prod_mov_p = exp(lambdaM./dotM);
                sum_mov_p  = lambdaM*difM/dotM;
            end
        end
        
        % Check if we are too close into a pole
        delta_x = x(:)-xold1(:) ;
        if 1e4 < abs( prod_fix_p * prod_mov_p )
            disp('abnormal obj values');
            x(:) = xold1(:) + delta_x/2;
            x(x<0) = 0 ;
            x(x>1) = 1 ;
            break ;
%             continue;
        end
        step = norm(delta_x) ;
        dist = dot(x(:)-x_star(:),x(:)-x_star(:))/(nelx*nely) ;
        
        %% Compute tunneling function
        if strcmp(tunneling,'classical') 
            f0val   = (c-fs)/(prod_fix_p*prod_mov_p);
            df0dx   = 1/prod_fix_p*(dc(:)-2*(c-fs)*(sum_fix_p+sum_mov_p));
            df0dxNP = prod_fix_p*(dc(:)-2*(c-fs)*(sum_fix_p));
        elseif strcmp(tunneling,'exponential') 
            f0val   = (c-fs)*prod_fix_p*prod_mov_p ;
            df0dx   = prod_fix_p*(dc(:)-2*(c-fs)*(sum_fix_p+sum_mov_p));
            df0dxNP = prod_fix_p*(dc(:)-2*(c-fs)*(sum_fix_p));
        elseif strcmp(tunneling,'gaussian') 
            f0val = c*(prod_fix_p*prod_mov_p);
            df0dx = (prod_fix_p*prod_mov_p) * ( dc(:) - 2*c*(sum_fix_p+sum_mov_p) );
%             df0dxNP = prod_fix_p*(dc(:)-2*(c-fs)*(sum_fix_p));
        end
        
        if lambdaM == 0 
            a = 0 ;
        else
            a = 1 ;
        end
        
        % Print to console and to file
        fprintf(['TS:%3i  It:%4i  Obj:%7.3f  Cons:%0.3f  KKT-norm.:%7.3f  RChange: %4.3f  T(x):%5.3f  Step:%7.3f  Dist:%0.4f  p:TS:%1i\n'], ...
            tunneling_step,loop,c,mean(xPhys(:)),kktnorm,rel_obj_change,f0val,step,dist,a);
%         fprintf(fileID,'%5i,%5i,%11.4f,%7.3f,%11.4f,%11.4f\n', ...
%             tunneling_step,loop,c,mean(xPhys(:)),f0val,norm(x(:)-x_star(:)));
        hist(counter,:) = [ tunneling_step,loop,c,mean(xPhys(:)),kktnorm,rel_obj_change,f0val,step,dist,a] ;

        if 0 < loop
            rel_obj_change = abs((f0val - old_obj)/old_obj);
        end
        
        %% Optimize
        if 0 < loop
            [~,kktnorm,~] = ...
                kktcheck(ncons,ndv,xmma,ymma,zmma,lam,xsi,eta,muu,zet,s, ...
                xmin,xmax,df0dx,fval,dfdx,a0,ai,ci,d);
        end
        
        [xmma,ymma,zmma,lam,xsi,eta,muu,zet,s,low,upp] = ...
            mmasub(ncons,ndv,loop,xval,xmin,xmax,xold1, ...
            xold2, f0val,df0dx,fval,dfdx,low,upp,a0,ai,ci,d);
        xnew     = reshape(xmma,nely,nelx);
        
        % Modify densities
        xPhys = filterDensities(xnew) ;
        % xProj = projectDensities(xPhys) ;
        xProj = xPhys;
        
        % Update MMA history
        xold2   = xold1(:);
        xold1   = x(:);
        x       = xnew;
        old_obj = f0val;
        
        % Plot densities
        colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        
        loop = loop + 1;
        counter = counter + 1 ;
    end
    
    % Finalization criterium
    if(fs-c>0)
        disp('Tunneling found root!');
        break;
    else
        xp = [xp,x(:)] ;
        disp('Tunneling failed, adding minimum');
    end
end
if counter < max_ts*maxoutit
   hist = hist(1:counter-1,:) ; 
end
end