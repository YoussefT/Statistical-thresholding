function [ s ] = constantrate_lrtstat( xvect)

xmeanjk = mean(mean(xvect));
xmeank  = mean(xvect,1);

s = 2*size(xvect,1)*sum(xmeank.*log(xmeank/xmeanjk));

end