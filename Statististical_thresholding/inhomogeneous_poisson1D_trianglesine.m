function [xvect] = inhomogeneous_poisson1D_trianglesine (lambda0,T,V,nu,A)

max = fminbnd(@(x) -trianglesine_fun(x,lambda0,T,V,nu,A),0.2*T,0.4*T);
lambdamax = trianglesine_fun(max,lambda0,T,V,nu,A);
xvecth = homogeneous_poisson1D( lambdamax,T);
n=length(xvecth);
lambdavect=trianglesine_fun(xvecth,lambda0,T,V,nu,A);
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

