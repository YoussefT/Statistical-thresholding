function [xvect] = inhomogeneous_poisson1D_blocks (lambda0)

T = 1;
max = fminbnd(@(x) -blocks_fun(x,lambda0),0.5,T);
lambdamax = blocks_fun(max,lambda0);
xvecth = homogeneous_poisson1D( lambdamax,T);
n=length(xvecth);
lambdavect=blocks_fun(xvecth,lambda0);
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