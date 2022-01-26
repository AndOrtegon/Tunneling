function x = randomDirection(x_star, epsilon)
    [nely, nelx] = size(x_star) ;
    ndv = nelx * nely ;

    rand_dir = 2*rand(1,ndv) - 1; %random number between -1 and 1
    rand_dir = reshape(rand_dir,nely,nelx);
    x = x_star + epsilon*rand_dir;

%     rand_dir = rand(1,ndv) ;
%     rand_dir(x_star > (0.5)) = -rand_dir(x_star > (0.5)) ;
%     rand_dir = rand_dir/norm(rand_dir);
%     rand_dir = reshape(rand_dir,nely,nelx);
%     x = x_star + epsilon*rand_dir;

    x(x>1) = 1;
    x(x<0) = 0;
end