close all
reset_random;

normF = @(x) norm(x,'fro');

%---------------------------------------------------------------------
% Prepare ground truth
n_imm = 8;  % # of images
noise = 1e-4; % gaussian noise on image points
density = .8; % holes in the adjacency and visibility

% 3D scene
d = 8; % distance camera - object center

% 3D points
load tribuna;
n_pts = length(vertices);
M =  vertices + [0;0;1];

% cameras
% internal parameters
width = 480; height=360;
K = par2K([width/2,height/2, -1.4  1 0]);

G_gt = cell(1,n_imm);
P_gt = cell(1,n_imm);
C_gt = cell(1,n_imm);
m = cell(1,n_imm);

for i = 1:n_imm
    C_gt{i} =  [5 * (rand(2,1) -.5);d+ 2*(rand-.5)]; % cop
    G_gt{i} = camera(  C_gt{i}, rand(3,1)-.5, [0;1;0]);
    P_gt{i}  = K * G_gt{i};
    G_gt{i} = [G_gt{i}; 0 0 0 1];
    
    % image points
    m{i} = htx(P_gt{i},M)  + noise*randn(2,n_pts); % noise
end

% some plots
scatter3(M(1,:), M(2,:),M(3,:),[],lines(n_pts),'filled');
hold on, wireframe(M', edges, 'k:'); hold on

for i = 1:n_imm
    plotcam(P_gt{i},.7);
end

disp(' ');
%-------------------------------------------------------------------------
% SMZ calibration

% planar points on a grid
[Xgrid,Ygrid] = meshgrid(-0.5:0.2:0.4);
M_grid = [Xgrid(:)';Ygrid(:)'; zeros(1,size(Xgrid(:),1))];
plot3(M_grid(1,:), M_grid(2,:), M_grid(3,:), '.r')
zlabel('Z')

H_est=cell(1,1);
m_grid =  cell(1,n_imm);
for i = 1:n_imm
    m_grid{i}  = htx(P_gt{i},M_grid)  + noise*randn(size(2,100)); % noise
    
    % H{i} = dlt(m{i},M);
    H_est{i}  = hom_lin(m_grid{i}, M_grid(1:2,:));
end

[P_est ,K_est] = calibSMZ(H_est);

fprintf('CalibSMZ reproj RMS error:\t\t %0.5g \n',...
    rmse(reproj_res_batch(P_est,M_grid, m_grid)) );

[P_out,M_out] = bundleadj(P_est,M_grid,m_grid,...
    'AdjustCommonIntrinsic','IntrinsicParameters',5, 'FixedPoints',size(M_grid,2));

fprintf('CalibSMZ reproj RMS error:\t\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,M_out,m_grid)) );


%--------------------------------------------------------------
% (Auto)Calibration from H_infty LVH

H12= K*eul(rand(1,3))*inv(K);
H13= K*eul(rand(1,3))*inv(K);
H23= K*eul(rand(1,3))*inv(K);

K_out  = calibLVH( {H12,H13,H23 });

fprintf('H_infty calibration %% error:\t\t %0.5g \n',100*abs((K_out(1,1)-K(1,1))/K(1,1)) );

%-------------------------------------------------------------------------
% Triangulation


% random visibility
while true
    % if density is too low this might cicle forever
    vis = rand(n_pts,n_imm) < density; % is logical
    figure, spy(vis),title('Visibility');ylabel('points');xlabel('images')
    if all(sum(vis,2) >= 3)
        break
    end
end

X_est=triang_lin_batch(P_gt, m, vis);
fprintf('Triangulation batch error:\t\t %0.5g \n', normF(M-X_est)/n_pts  );

X_est = triang_nonlin_batch(X_est,P_gt, m, vis);
fprintf('Triangulation non_lin batch error:\t %0.5g \n', normF(M-X_est)/n_pts  );

%-------------------------------------------------------------------------
% Projective reconstruction

[P_proj,M_proj] = prec(m);

fprintf('Projective Recon reproj RMS error:\t %0.5g \n',...
    rmse(reproj_res_batch(P_proj,M_proj, m)));

%---------------------------------------------------------------------
% Bundle Adjustment

kappa = num2cell(.2*ones(1,n_imm),1);
M_in = X_est;
% Initial rec
for i = 1:n_imm
    P_in{i} = P_gt{i};
    m{i} = rdx(kappa{i},m{i},K);  % add radial
end

fprintf('Bundle Adjustment RMS error (before):\t %0.5g \n',...
    rmse(reproj_res_batch(P_in, M_in, m, 'Visibility', vis)) );

%[P_out,M_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'FixedIntrinsic');
[P_out,M_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'FixedIntrinsic', 'DistortionCoefficients', kappa);
kappa_out = kappa;
fprintf('Bundle Adjustment RMS error (after):\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,M_out, m, 'Visibility',vis,'DistortionCoefficients', kappa_out) ));


% [P_out,M_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'AdjustCommonIntrinsic','IntrinsicParameters',5);
[P_out,M_out, kappa_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'AdjustCommonIntrinsic','IntrinsicParameters',5,'DistortionCoefficients', num2cell(zeros(1,n_imm),1));
fprintf('Bundle Adjustment RMS error (after):\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,M_out, m, 'Visibility',vis,'DistortionCoefficients', kappa_out) ));


% [P_out,M_out,kappa_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'AdjustSeparateIntrinsic');
[P_out,M_out,kappa_out] = bundleadj(P_in,M_in,m,'Visibility',vis,'AdjustSeparateIntrinsic','DistortionCoefficients', num2cell(zeros(1,n_imm),1));
fprintf('Bundle Adjustment RMS error (after):\t %0.5g \n', ...
    rmse(reproj_res_batch(P_out,M_out, m, 'Visibility',vis,'DistortionCoefficients', kappa_out) ));

disp(' ');
%---------------------------------------------------------------------
%% Sinchronization in SE(3)

% Compute relative orientations
X = cell2mat(G_gt(:));
Y = cell2mat(cellfun(@inv,G_gt(:)','uni',0));
Z = X * Y;

% random adjacency matrix
A = rand(n_imm) < density;
A = triu(A,1) + triu(A,1)';
figure, spy(A), xlabel(''); title('Adjacency SE(3)');

% make holes in Z
Z = Z.*kron(A,ones(4));

% is connected?
if any(any((eye(n_imm)+A)^n_imm==0))
    error('graph is not connected, synchronization is impossible')
end

%---------------------------------------------------------------------
%% Rotation synch

% squeeze out rotations from Z
Z_rot = Z;
Z_rot(4:4:end,:)=[]; Z_rot(:,4:4:end)=[];
R = rotation_synch(Z_rot,A);

%---------------------------------------------------------------------
%% Translation synch

%  gather bearings from Z and compute baselines using R
U = extract_bearings(Z,R) ;
B = adj2inc(A);
C = translation_synch(U,B);

% compute error in SE(3) (rotation+translation)
err =0;
for i=1:n_imm
    Gi=[R{i},-R{i}*C{i}; 0 0 0 1];
    err = err + normF(Gi*G_gt{1} - G_gt{i}) ;
end

fprintf('Synchronization error (SE3):\t %0.5g \n',err/n_imm);

%---------------------------------------------------------------------
%% Localization from bearings

% is paralell-rigid?
if ~ParallelRigidityTest(A,3)
    error('graph is not parallel rigid, localization is impossible')
end

% forget the module of translation
for k=1:size(U,2)
    U(:,k) = U(:,k)/norm(U(:,k));
end

C = cop_from_bearings(U,B);

% [R,t,s] = opa(reshape(cell2mat(C_gt'),3,[]),reshape(cell2mat(C'),3,[]));
R = G_gt{1}(1:3,1:3)';
s = norm(C_gt{2}-C_gt{1})/norm(R*C{2});
t = C_gt{1}/s;

% compute error
err =0;
for i = 1:n_imm
    err = err + norm(C_gt{i} - s*(R*C{i}+t));
end

fprintf('Localization error:\t\t %0.5g \n',err/n_imm);


%--------------------------------------------------------------
%% Autocalibration

% obtain F from Z
F=cell(1,1);
for j=1:size(A,1)
    for i=j+1:size(A,2)
        if A(i,j) > 0
            F{i,j} = inv(K)'*skew(Z(4*i-3:4*i-1,4*j))*Z(4*i-3:4*i-1,4*j-3:4*j-1)*inv(K);
        end
    end
end

K0 = K + 20*randn(3,3); K0(3,3) = 1;
K_out = autocal(F,K0)

fprintf('Autocalibration %% error:\t %0.5g \n',100*abs(K_out(1,1)-K(1,1))/abs(K(1,1))) ;

%--------------------------------------------------------------
%% Homography synchronization

Hgt=cell(1,n_imm);
for i=1:n_imm
    Hgt{i} = randn(3);
    Hgt{i} = Hgt{i}./nthroot(det(Hgt{i}),3);
end

X = cell2mat(Hgt(:));
Y= cell2mat(cellfun(@inv,Hgt(:)','uni',0));
Z = X * Y;

% random adjacency matrix
A = rand(n_imm) < density;
A = triu(A,1) + triu(A,1)';

% make holes
Z = Z.*kron(A,ones(3));
figure, spy(A), xlabel(''); title('Adjacency H');

% is connected?
if any(any((eye(n_imm)+A)^n_imm==0))
    error('graph is not connected, synchronization is impossible')
end

H = hom_synch(Z,A);

err =0;
for i=1:n_imm
    err = err + normF(H{i}*Hgt{1} - Hgt{i}) ;
end

fprintf('Synchronization error (H):\t %0.5g \n',err/n_imm);





