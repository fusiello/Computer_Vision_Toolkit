close all
reset_random;

normF = @(x) norm(x,'fro');

%---------------------------------------------------------------------
% Prepare ground truth
noise = 1e-4;

% 3D scene
dist = 3; % distance camera - object center 

% 3D points
load tribuna;
n = length(vertices);
X = [0;0;dist] + vertices;
% X =  [0;0;dist] + 2*(rand(3, n) - .5);

% cameras
% internal parameters
width = 480; height=360;
K = par2K([width/3,height/3, -1.4  1 0]);
%P1 = K*[eye(3), [ -1;0;0]];
%P2 = K*[eye(3), [ 1;0;0]];
P1= K*camera([ .9;0;0],[-.05; .05; dist], [.95; 1; 0]); %left
P2= K*camera([-.9;0;0],[.05; -.05; dist], [-.97;1; 0]); %right

% image points
x1 = htx(P1,X)+ noise*randn(2,n);
x2 = htx(P2,X)+ noise*randn(2,n);

% some plots
figure,
scatter3(X(1,:), X(2,:),X(3,:),[],lines(size(x1,2)),'filled');
hold on, wireframe(X', edges, 'k:');
plotcam(P1,.5,'r'); plotcam(P2,.5,'b');
xlabel('X'), ylabel('Y'), zlabel('Z')
view(180, -90 )

figure;
subplot(1,2,1); scatter(x1(1,:),x1(2,:),[],lines(size(x1,2)),'filled');
hold on, wireframe(x1', edges, 'k:');
title('Image 1'), axis equal,ylim([0 height]); xlim([0 width]);
axis ij, set(gca,'XAxisLocation','top')

subplot(1,2,2); scatter(x2(1,:), x2(2,:),[],lines(size(x1,2)),'filled');
hold on, wireframe(x2', edges, 'k:');
title('Image 2'), axis equal, ylim([0 height]); xlim([ 0 width]),
axis ij,  set(gca,'XAxisLocation','top')

%-------------------------------------------------------------------------
% Rectification
[T1,T2,Pn1,Pn2] = rectifyP(P1,P2);
xr1 = htx(T1,x1); xr2 = htx(T2,x2);

% some plots
figure;
subplot(1,2,1); scatter(xr1(1,:), xr1(2,:),[],lines(size(x1,2)),'filled');
hold on, wireframe(xr1', edges, 'k:');
title('Image 1 rectified'), axis equal, ylim([0 height]); xlim([0 width]);
axis ij, set(gca,'XAxisLocation','top'),
set(gca, 'YMinorGrid', 'on', 'XGrid', 'off', 'MinorGridLineStyle', '-')
yt = yticks;

subplot(1,2,2); scatter(xr2(1,:), xr2(2,:),[],lines(size(x1,2)),'filled');
hold on, wireframe(xr2', edges, 'k:');
title('Image 2 rectified'),  axis equal,ylim([0 height]); xlim([ 0 width]),
axis ij,  set(gca,'XAxisLocation','top')
set(gca, 'YMinorGrid', 'on', 'XGrid', 'off', 'MinorGridLineStyle', '-')
yticks(yt);

figure(1), plotcam(Pn1,.5,'k'); plotcam(Pn2,.5,'k');

% model from disparity
X_model = triang_disp( xr1(1,:),xr1(2,:),xr2(1,:)-xr1(1,:), Pn1,Pn2 );

% align to GT
[R,t,s] = opa(X,X_model);
X_obj = s*(R*X_model + t*ones(1,size(X,2)));
fprintf('Triangulation f/ disparity err:\t %0.5g \n',  rmse(X(:)-X_obj(:)));

figure
plot3(X(1,:), X(2,:), X(3,:), 'or'); hold on
plot3(X_obj(1,:), X_obj(2,:), X_obj(3,:),'+b');
title('Triangulation from disparity')

disp(' ');

%--------------------------------------------------------------------------
% radial distortion
%
k_1 = .2;
[x_dist,J] = rdx(k_1,x1,K);

fprintf('Radial sanity test:\t\t %0.5g \n', norm(irdx(k_1,x_dist,K)- x1) );

figure
plot (x1(1,:),x1(2,:),'k.'); hold on
plot (x_dist(1,:),x_dist(2,:),'r.');

% Warning: this procedure is rather unstable, it can converge to wrong
% result depending on the distribution of points in the image
k = calib_radial(x_dist, X, 1);

fprintf('Radial calibation error:\t %0.5g \n', abs(k_1-k));
x_u = irdx(k,x_dist,K);
plot (x_u(1,:),x_u(2,:),'go');
legend('Ideal','Distorted','Corrected')

disp(' ');

%--------------------------------------------------------------------------
% resection

P_est = resect_lin(x1, X);
x_est = htx(P_est,X);  % project with estimated camera
fprintf('Resection ___lin reproj error:\t %0.5g \n', rmse(x1(:)-x_est(:)));

P_est = resect_nonlin(P_est, x1, X);
x_est = htx(P_est,X);  % project with estimated camera
fprintf('Resection nonlin reproj error:\t %0.5g \n', rmse(x1(:)-x_est(:)));

disp(' ');
%--------------------------------------------------------------------------
% Triangulation

X1 = triang_lin({P1, P2}, {x1(:,1) ,x2(:,1)});
fprintf('Triangulation 1pt ___lin error:\t %0.5g \n', normF(X(:,1)-X1) );

X1_nl = triang_nonlin(X1, {P1, P2}, {x1(:,1) ,x2(:,1)});
fprintf('Triangulation 1pt nonlin error:\t %0.5g \n', normF(X(:,1)-X1_nl) );

X_est = triang_lin_batch({P1, P2}, {x1,x2});
fprintf('Triangulation batch error:\t %0.5g \n', rmse(X(:)-X_est(:)));

disp(' ');
%-------------------------------------------------------------------------
% separate exterior orientation

[R1,t1] = exterior_lin(x1,X,K);
fprintf('Exterior ___lin SE3 error:\t %0.5g \n',normF([R1,t1] -  K\P1));
%
[R1,t1] = exterior_iter(x1,X,K);
fprintf('Exterior __iter SE3 error:\t %0.5g \n',normF([R1,t1] -  K\P1));

% refine 1
[R1,t1] = exterior_nonlin(R1, t1 , x1 , X,  K);
fprintf('Exterior nonlin SE3 error:\t %0.5g \n',normF([R1,t1] -  K\P1));

[R2,t2] = exterior_lin(x2,X,K);
fprintf('Exterior ___lin SE3 error:\t %0.5g \n', normF([R2,t2] -  K\P2));

% should refine 2 as well...

X_obj = triang_lin_batch({K*[R1,t1], K*[R2,t2]}, {x1, x2});
fprintf('Separate exterior triang error:\t %0.5g \n', rmse(X(:)-X_obj(:)));

figure
plot3(X(1,:), X(2,:), X(3,:), 'or'); hold on
plot3(X_obj(1,:), X_obj(2,:), X_obj(3,:),'+b');
title('Separate exterior o.')

%--------------------------------------------------------------------------
[P_out,X_out] = bundleadj({K*[R1,t1], K*[R2,t2]},X_obj,{x1, x2});
fprintf('Bundle adjustment RMS error:\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,X_out, {x1, x2})));

disp(' ');
%--------------------------------------------------------------------------
% relative orientation + absolute orientation
G21 = [K\P2; 0 0 0 1]  * inv([ K\P1; 0 0 0 1]);

[R21,t21] = relative_lin(x2, x1, K, K);
fprintf('Relative ___lin SO3 error:\t %0.5g \n', normF(R21 - G21(1:3,1:3)));

[R21,t21] = relative_nonlin(R21,t21 ,x2, x1, K, K);
fprintf('Relative nonlin SO3 errror:\t %0.5g \n',  normF(R21 - G21(1:3,1:3) ));

X_model = triang_lin_batch({K*[eye(3),zeros(3,1)], K*[R21,t21]}, {x1,x2});

% align to GT - assume the first 6 points are GCP
[R,t,s] = opa(X(:,1:6),X_model(:,1:6));
X_obj = s*(R*X_model + t*ones(1,size(X,2)));
fprintf('Relative triang error:\t\t %0.5g \n',  rmse(X(:)-X_obj(:)));

figure
plot3(X(1,:), X(2,:), X(3,:), 'or'); hold on
plot3(X_obj(1,:), X_obj(2,:), X_obj(3,:),'+b');
title('Relativo o.')
%--------------------------------------------------------------------------

[P_out,X_out] = bundleadj({K*[R1,t1], K*[R2,t2]},X_obj,{x1, x2});
fprintf('Bundle adjustment RMS error:\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,X_out, {x1, x2})));

disp(' ');
%-------------------------------------------------------------------------
% Essential

E = skew(G21(1:3,4))* G21(1:3,1:3);

E_est = essential_lin(x2,x1,K,K);
fprintf('Essential Smps error:\t\t %0.5g \n', rmse(sampson_fund(inv(K)'*E_est*inv(K),x1,x2)));

disp(' ');
%-------------------------------------------------------------------------
% Fundamental

F = fund(P1,P2);

F_est = fund_lin(x2,x1);
fprintf('Fundamental Smps error:\t\t %0.5g \n', rmse(sampson_fund(F_est,x1,x2)));

F_out = fund_nonlin(F_est, x2, x1);
fprintf('Fundamental nonlin Smps error:\t %0.5g \n', rmse(sampson_fund(F_out,x1,x2)));

out=zeros(size(x1)); out(:,1) = 100;
[F_msac, in]  = fund_rob(x2,x1+out,'MSAC',1);
fprintf('Fundamental MSAC Smps error:\t %0.5g \n', rmse(sampson_fund(F_msac,x1(:,in),x2(:,in))));

[F_lms,in]  = fund_rob(x2,x1+out,'LMS');
fprintf('Fundamental LMS Smps error:\t %0.5g \n', rmse(sampson_fund(F_lms,x1(:,in),x2(:,in))));

[F_irls,in]  = fund_rob(x2,x1+out,'IRLS');
fprintf('Fundamental IRLS Smps error:\t %0.5g \n', rmse(sampson_fund(F_irls,x1(:,in),x2(:,in))));


disp(' ');
%-------------------------------------------------------------------------
% Homography

[Xgrid,Ygrid] = meshgrid(-.5:0.2:.5);
M_grid = [Xgrid(:)';Ygrid(:)'; 2*ones(1,size(Xgrid(:),1))];
figure(1), hold on, plot3(M_grid(1,:), M_grid(2,:), M_grid(3,:), '.')

x1_grid  = htx(P1,M_grid) + noise*randn(size(M_grid)- [1 0]);
x2_grid  = htx(P2,M_grid) + noise*randn(size(M_grid)- [1 0]);
n_grid = size(x1_grid,2);

H_est  = hom_lin(x2_grid,x1_grid);
fprintf('Homography ___lin Smps error:\t %0.5g \n', rmse(sampson_hom(H_est,x1_grid,x2_grid)));

H_out = hom_nonlin(H_est, x2_grid, x1_grid);
fprintf('Homography nonlin Smps error:\t %0.5g \n', rmse(sampson_hom(H_out,x1_grid,x2_grid)));

%%
out=zeros(size(x1_grid)); out(:,1) = 10;
[H_msac, in]  = hom_rob(x2_grid+out,x1_grid,'MSAC',1);
fprintf('Homography MSAC Smps error:\t %0.5g \n', rmse(sampson_hom(H_msac,x1_grid(:,in),x2_grid(:,in))));

[H_msac, in]  = hom_rob(x2_grid+out,x1_grid,'LMS');
fprintf('Homography LMS Smps error:\t %0.5g \n', rmse(sampson_hom(H_msac,x1_grid(:,in),x2_grid(:,in))));

[H_msac, in]  = hom_rob(x2_grid+out,x1_grid,'IRLS');
fprintf('Homography IRLS Smps error:\t %0.5g \n', rmse(sampson_hom(H_msac,x1_grid(:,in),x2_grid(:,in))));

%%
H_est  = hom_lin(x2_grid,M_grid(1:2,:));

fdh(H_est,K(1,3),K(2,3));

%%

disp(' ');
%-------------------------------------------------------------------------
% Uncalibrated

T = rand(4,4); % projective transform
P1p = P1 * T;
P2p = P2 * T;

[P1o,P2o,T] = eucl_upgrade(P1p, P2p, K, K);
P1_est = P1p * T;
P2_est = P2p * T;

% norm(P2est(:,1:3) - P2o(:,1:3))

G21_est = [K\P2_est; 0 0 0 1]  * inv([ K\P1_est; 0 0 0 1]);
fprintf('Euc upgrade SO3 error:\t\t %0.5g \n', normF(G21_est(1:3,1:3) - G21(1:3,1:3)));

