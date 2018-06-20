function [ test, pvalue, R ] = innovation_test(points,T,L,alpha)

if not(iscell(points))
    point_cell{1} = points;
else
    point_cell = reshape(points, max(size(points)),1);
end

M = size(point_cell,1);
c = zeros (M,2^(L+1));

for m=1:M
    c(m,:) = haar_counts(point_cell{m},T,L+1);
end
R = pairconstrate_lrtstat(c,2^L);
test = R > chi2inv(1-alpha,2^(L)) ;
pvalue = chi2cdf(R,2^(L),'upper');


end