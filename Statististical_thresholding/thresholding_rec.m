function [ rate_est ] = thresholding_rec( points,t,T,j0,J,p )

%%% Estimate initial projection on V_j0
rate_est = thresholding_linear(points,t,T,j0);


betas_rec = zeros(J-j0+1,2^J);
M = size(points,1);
counts = cell(M,J-j0+1);
for j=j0:J
    for m = 1:M
        counts{m,j-j0+1} = haar_counts(points{m},T,j+1);
    end
end

for j=j0:J
    counts2mat = cell2mat(counts(:,j-j0+1));
    beta_j = zeros(1,2^j);
    idx_list = [];
    for k = 1:2^j
        beta_j(k) = (2^(j/2))*(mean(counts2mat(:,2*k-1),1)-mean(counts2mat(:,2*k),1))/sqrt(T);         
    end
    P = 2^j;
    R = pairconstrate_lrtstat(counts2mat,P);
    test = R > chi2inv(p,2^j);
    while test 
       idx = find (abs(beta_j) == max(abs(beta_j)),1);
       idx_list = [idx_list, 2*idx-1,2*idx];
       betas_rec(j-j0+1,idx) = beta_j(idx);
       beta_j(idx) = 0;
       counts_temp = counts2mat;
       counts_temp(:,idx_list) = [];
       P = P-1;
       R = pairconstrate_lrtstat(counts_temp,P);
       test = R > chi2inv(p,P);
    end
    S = repmat(betas_rec(j-j0+1,1:2^j)',1,length(t));
    H = squeeze(haar_base(t,0,T,j,j));
    if j == 0
        H = H';
    end
    w = S.*H;
    w = sum(w,1);
    rate_est = rate_est + w;
end

end

