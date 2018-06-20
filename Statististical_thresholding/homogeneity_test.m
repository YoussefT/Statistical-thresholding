function [ test, pvalue, R ] = homogeneity_test(points,T,J,alpha)

%Performs the likelihood ratio test for J-th level homogeneity on point process data.

if not(iscell(points))
    point_cell{1} = points;
else
    point_cell = reshape(points, max(size(points)),1);
end

M = size(point_cell,1);
c = zeros (M,2^(J));

for m=1:M
    c(m,:) = haar_counts(point_cell{m},T,J);
end
R = constantrate_lrtstat(c);
test = R > chi2inv(1-alpha,2^(J)-1) ;
pvalue = chi2cdf(R,2^(J)-1,'upper');


end
