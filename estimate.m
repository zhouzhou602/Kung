clear;

global delta;

% import household data
household_data = importdata('household_data.csv');
loginc = household_data.data(:,2);
college = household_data.data(:,3); 
z = [(loginc-mean(loginc))./std(loginc), college];  % normalize loginc to have mean zero and std 1
znames = household_data.textdata(2:3);
jchoice = household_data.data(:,4);
N = length(jchoice);
zk = length(znames);

% import neighborhood data
neighborhood_data = importdata('neighborhood_data.csv');
schqual = neighborhood_data.data(:,2);
price = neighborhood_data.data(:,3);
sharej = neighborhood_data.data(:,4);
iv = neighborhood_data.data(:,5);
x = [schqual, price];
xnames = neighborhood_data.textdata(2:3);
J = length(price);
xk = length(xnames);

% set the log likelihood function
obj = @(params)loglike(params, z, x, sharej, jchoice);

% initialize starting guess
B0 = zeros(zk, xk);
params0 = reshape(B0, zk, xk);
tic
obj0 = obj(params0);
toc

% estimate B
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',10000,'gradobj','on');
tic
paramshat = fminunc(obj, params0, options);
toc
B = reshape(paramshat, zk, xk);

% calculate standard errors for B
tic
[objhat, gradobj, avar] = obj(paramshat);
toc
Bse = reshape(diag(sqrt(avar)), zk, xk);

% estimate A
Y = delta;
X = [ones(J,1), x(:,1), x(:,2)];
Z = [ones(J,1), x(:,1), iv];
A = (Z'*X)\(Z'*Y);

Yhat = X*A;
bfs = (Z'*Z)\(Z'*x(:,2));
Xhat = [ones(J,1), x(:,1), Z*bfs];
rmse = sqrt(mean((Y-Yhat).^2));
Ase = rmse*sqrt(diag(inv(Xhat'*Xhat)));

est_all = [A(2:3)'; B];
se_all = [Ase(2:3)'; Bse];
mytab = esttab(est_all, se_all, 'estimates', [{'const'}, znames], xnames, {''}, {''});
showtable(mytab)

save('results_latest.mat');









