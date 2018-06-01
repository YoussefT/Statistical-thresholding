function [ rate_est ] = thresholding_linear( points,t,T,J )

M = size(points,1);
counts = zeros(M,2^J);
for m = 1:M
    counts(m,:) = haar_counts(points{m},T,J);
end
coeff = 2^(J/2)*mean(counts,1)/sqrt(T);
S = repmat(coeff',1,length(t));
H = reshape(squeeze(haar_base_father(t,0,T,J)),2^J,length(t));
w = S.*H;
rate_est = sum(w,1);


end

