
disp('Testing aux funxtions...')

a = (rand(3,1)); b = (rand(3,1));

% sanity check for eul/ieul
R = eul(a);
norm(R-eul(ieul(R)))

% sanity check for rod/irod
[u,theta] = irod(R); 
norm(R-rod(u,theta))

% sanity check for quat/iquat
norm(R-quat(iquat(R)))

% sanity check for skew/iskew
S = skew(a);
norm(S-skew(iskew(S)))

% test simpleGN (result should be [1,1])
vfun = @(x)[10*(x(2) - x(1)^2),1 - x(1)]';
jac = @(x)[-20*x(1),10; -1,0];        
vfunjac = @(x)deal(vfun(x),jac(x));
x0 = [0,2]';
x = simpleGN(vfunjac,x0)
x = lsq_nonlin(vfunjac,x0)


