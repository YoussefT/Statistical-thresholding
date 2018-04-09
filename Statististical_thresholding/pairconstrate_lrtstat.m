function [ s ] = pairconstrate_lrtstat( xvect,P)

xmeankm = mean(xvect,1);
xmeank  = zeros(size(xvect,2)/2,1);
for k = 1:P
    xmeank(k) = mean(mean(xvect(:,2*k-1:2*k)));
end
xmeank = kron(xmeank',ones(1,2));
logx = xmeankm.*log(xmeankm./xmeank);
logx(isnan(logx)) = 0;
s = 2*size(xvect,1)*sum(logx);

end