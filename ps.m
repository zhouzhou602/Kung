function loglike
B = [1 1;1 1]

%% run interation.m for deltafin
deltafin=iteration(sharej,B,jchoice,z,x); %%This is our list of deltafin
deltamatfin = repmat(deltafin,1,N)';
ufin = z*B*x' + deltamatfin; %%Construct a N by J matrix for u
exp_ufin = exp(ufin); %%Take the exponential function 
sum_exp_ufin_j = sum(exp_ufin,2); %%This is a N by 1 matrix
sum_exp_ufin = repmat(sum_exp_ufin_j,1,J); 
pfin = exp_u./sum_exp_u;

%%Then we would like to construct p* where the individual made their actual
%%choice of j.
choice = zeros(N,J)
jchoice = household_data.data(:,4);
for i = 1:N,
    choice(i,jchoice(i,1))=1;
end
p_real = pfin.*choice;
sum_rows_p = sum(p_real,2);
log_p_real = log(sum_rows_p);
LL = -sum(log_p_real,1); %%we got the likelihood