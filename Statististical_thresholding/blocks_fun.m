function [ flambda ] = blocks_fun(x,lambda0)

Iblocks = 1.551;

xj = [0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81];
hj = [4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2];

X = repmat(x,11,1);
Xj = repmat(xj',1,length(x));
Hj = repmat(hj',1,length(x));

flambda(:,1) = sum(Hj.*(1+sign(X-Xj))/2,1)/Iblocks;
flambda = 0.25*lambda0*flambda + 1.75*lambda0;

end
