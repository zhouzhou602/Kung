
function delta_prime = iteration(params,z,x,sharej,jchoice)
N = length(jchoice);
J=length(sharej);
delta1=sharej;
deltamat = repmat(delta1,1,N)';
u = z*params*x' + deltamat; %%Construct a N by J matrix for u
exp_u = exp(u); %%Take the exponential function 
sum_exp_u_j = sum(exp_u,2); %%This is a N by 1 matrix
sum_exp_u = repmat(sum_exp_u_j,1,J); 
frac = exp_u./sum_exp_u;
sum_frac = sum(frac,1); %%sum up over N 
share = (1/N)*sum_frac'; %%to get a J by 1 matrix for share

stop = 1e-09;

t = 1;
old_share = sharej;
while t < 2000
    delta_prime = delta1 + log(old_share) - log(share);
    diff = max(abs(delta_prime-delta1));
       if diff < stop
           break
       else
           delta1= delta_prime;
           t=t +1;
       end
end

 deltafin = delta_prime;

 
 