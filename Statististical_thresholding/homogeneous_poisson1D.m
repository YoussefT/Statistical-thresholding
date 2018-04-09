function [xvect] = homogeneous_poisson1D (lambda,X)


n = poissrnd(lambda*X);
if (n>0)
    xvect=rand(1,n)*X;
else
    xvect  = [];
end
xvect = sort(xvect);
% plot(sort(xvect),zeros(n,1),'.')
% xlabel('X')

end

