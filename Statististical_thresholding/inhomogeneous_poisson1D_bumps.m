function [xvect] = inhomogeneous_poisson1D_bumps (lambda0)

T = 1;
% max = fminbnd(@(x) -bumps_fun(x,lambda0),0.75*T,0.85*T);
% lambdamax = bumps_fun(max,lambda0);
lambdamax = bumps_fun(0.78,lambda0);
xvecth = homogeneous_poisson1D( lambdamax,T);
n=length(xvecth);
lambdavect=bumps_fun(xvecth,lambda0);
p=lambdavect/lambdamax;
i=1;
xvect = zeros(n,1);    

while (i<=n)
u=rand;
if (u<p(i))
xvect(i)=xvecth(i);
end
i=i+1;
end
xvect=xvect(xvect>0);

end