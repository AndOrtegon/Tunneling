function [x_star,f_star,hist] = optimization_phase(x,fileID,prog_penal)
    global nelx nely volfrac penal  ml penal_type plotting
    hist=0;
    % Continuation parameters
%     prog_penal = true;
    
    step_p = 0.02;
    target_p = penal;
    
    if strcmp(penal_type, 'SIMP')
        init_p = 1.0;
    elseif strcmp(penal_type, 'RAMP')
        init_p = 0.0;
    end    
    

    %% Modify densities
    xPhys = filterDensities(x) ;
    xProj = projectDensities(xPhys);

    %% INITIALIZE MMA PARAMETERS
    change = 1;
    outeriter = 0;
    maxoutit  = 500;

%     hist = zeros(maxoutit-1,nely*nelx) ;
    
    kkttol  = 1e-4;
    obj_tol = 1e-5;

%     TOLF = 1e-5;
     TOLF = 1e-7;
    ndv = size(x(:),1);
    ncons = 1;

    a0 = 1;
    ai = zeros(ncons,1);
    ci = 1000*ones(ncons,1);
    d = 0;

    low = ones(ndv,1);
    upp = ones(ndv,1);

    xold1=x(:);
    xold2=x(:);

    %% START ITERATION
    %%%% The iterations start:
    kktnorm = kkttol+10;
    loop = 0;
    old_obj = 0;
    rel_obj_change = 10*obj_tol;

    while loop < maxoutit
        %% CONTINUATION
        if prog_penal
            if loop == 0
                penal = init_p;
            else
                penal = min(penal + step_p, target_p);
            end
        end

        %% Optimization
        [c, vol, dc, dv] = analyze(xPhys,xProj) ;
        
        xval = x(:);
        xmin = max(0,xval-ml);
        xmax = min(1,xval+ml);

        f0val = c;
        df0dx = dc(:);
        fval  = vol/(volfrac*nely*nelx) - 1;
        dfdx  = dv(:)'/(volfrac*nely*nelx);

        fprintf('Optm    It:%4i  Obj:%7.3f  Cons: %4.3f  KKT-norm.:%7.3f  RChange: %4.3f\n', ...
                loop,c,mean(xPhys(:)),kktnorm,rel_obj_change);
                
        if 0 < loop
            rel_obj_change = abs((f0val - old_obj)/old_obj);
            
            if(abs(f0val-old_obj)<=TOLF*(1+f0val) && penal >= target_p)
                disp('TOLF requirement satisfy')
                break
            end
            
%             hist(loop+1,:) = x(:) - xold1(:) ;
            
            %       if rel_obj_change < obj_tol
            %           disp('Objective convergence tolerance satisfied.');
            %       end

            %     if(change < 0.01 && penal >= target_p)
            %         disp('DV change requirement satisfied.')
            %     end
        end

        % The residual vector of the KKT conditions is calculated:
        if 0 < loop
            [~,kktnorm,~] = ...
                kktcheck(ncons,ndv,xmma,ymma,zmma,lam,xsi,eta,muu,zet,s, ...
                xmin,xmax,df0dx,fval,dfdx,a0,ai,ci,d);
        end

        %%%% The MMA subproblem is solved at the point PARAM_VALUE
        [xmma,ymma,zmma,lam,xsi,eta,muu,zet,s,low,upp] = ...
            mmasub(ncons,ndv,loop,xval,xmin,xmax,xold1, ...
            xold2, f0val,df0dx,fval,dfdx,low,upp,a0,ai,ci,d);
        
        xnew     = reshape(xmma,nely,nelx);

        % Filtering and projection
        xPhys = filterDensities(xnew) ;
        xProj = projectDensities(xPhys);

        % Update MMA
        xold2 = xold1(:);
        xold1 = x(:);
        x     = xnew;

        old_obj = f0val;
        f_star = c;
        
        % % PRINT RESULTS
        %   fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        %       mean(xPhys(:)),change);
        fprintf(fileID,'%5i,%11.4f,%7.3f\n',loop,c,mean(xPhys(:)));
        
        % PLOT DENSITIES
        if plotting
            figure(1) 
            title('Optimization') 
            subplot(1,1,1)
            colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        end
        
        loop = loop + 1;
    end
    x_star = x;

    % Compute and report GRF
    sizex = numel(x_star);
    x_star_vec = reshape(x_star, sizex,1);

    % This assumes elements are of size 1 x 1, which is the case in this
    % implementation
    GRF = (4/sizex)*dot(x_star_vec,(1-x_star_vec));
    disp('GRF fraction for optimal design:');
    disp(GRF);

end