
disp('Testing aux funxtions...')

R = eul(rand(1,3));

norm(R-eul(ieul(R)))

[u,theta] = irod(R); norm(R-rod(u,theta))

norm(R-quat(iquat(R)))

S = skew(rand(1,3));

norm(S-skew(iskew(S)))


a = (rand(1,3));
b = (rand(1,3));
mu = rand;
lambda= rand;
c = mu*a + lambda*b;

[ mu_out, lambda_out ] = icomb( a, b, c );

mu - mu_out
lambda - lambda_out