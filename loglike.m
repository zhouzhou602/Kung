% run interation.m for deltafin
function [loglike,gradobj,avar] = loglike(params,z,x,sharej,jchoice) %%This is our list of deltafin

N = length(jchoice);
J = length(sharej);


deltafin = iteration(params,z,x,sharej,jchoice);
deltamatfin = repmat(deltafin,1,N)';
ufin = z*params*x' + deltamatfin; %%Construct a N by J matrix for u
exp_ufin = exp(ufin); %%Take the exponential function, still N by J
sum_exp_ufin_j = sum(exp_ufin,2); %%This is a N by 1 matrix
sum_exp_ufin = repmat(sum_exp_ufin_j,1,J); 
pfin = exp_ufin./sum_exp_ufin; %%N by J

%%Then we would like to construct p* where the individual made their actual
%%choice of j.
choice = zeros(N,J);
for i = 1:N 
    choice(i,jchoice(i,1))=1;
end
p_real = pfin.*choice; %%N by J
sum_rows_p = sum(p_real,2);
log_p_real = log(sum_rows_p);
loglike = -sum(log_p_real,1);%%we got the likelihood
z1x1 = z(:,1)*x(:,1)'; %%N by J
z1x2 = z(:,1)*x(:,2)';
z2x1 = z(:,2)*x(:,1)';
z2x2 = z(:,2)*x(:,2)';
kung11 = sum(z1x1.*choice,2)-sum(exp_ufin.*z1x1,2)./sum_exp_ufin_j;
sumkung11 = -sum(kung11,1);
kung12 = sum(z1x2.*choice,2)-sum(exp_ufin.*z1x2,2)./sum_exp_ufin_j;
sumkung12 = -sum(kung12,1);
kung21 = sum(z2x1.*choice,2)-sum(exp_ufin.*z2x1,2)./sum_exp_ufin_j;
sumkung21 = -sum(kung21,1);
kung22 = sum(z2x2.*choice,2)-sum(exp_ufin.*z2x2,2)./sum_exp_ufin_j;
sumkung22 = -sum(kung22,1);
gradobj = [sumkung11 sumkung21 sumkung12 sumkung22]'; %%4 by 1
score = [kung11 kung21 kung12 kung22]'; %%4 by 1
score_sq = score*score';
info_hat = (1/N)*score_sq; %%4 by 4
avar = (1/N)*inv(info_hat);