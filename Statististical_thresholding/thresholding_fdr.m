function [rate_est] = thresholding_fdr( points,t,T,j0,J,q )

%%% Estimate initial projection on V_j0
rate_est = thresholding_linear(points,t,T,j0);


pvalue = zeros(2^(J+1)-1,3);
m = 2^(J+1) - 2^(j0);
idx = 1:1:m ;
qm = q/sum(idx.^-1);
M = size(points,1);
counts = cell(M,J-j0+1);
for j=j0:J
    for m = 1:M
        counts{m,j-j0+1} = haar_counts(points{m},T,j+1);
    end
end

for j=j0:J
    counts2mat = cell2mat(counts(:,j-j0+1));
    for k = 0:2^j-1
        R = pairconstrate_lrtstat(counts2mat(:,2*k+1:2*k+2),1);
        pvalue(2^j+k,1) = chi2cdf(R,1,'upper');
        pvalue(2^j+k,2) = j;
        pvalue(2^j+k,3) = k;
    end        
end

pvalue = pvalue(2^(j0):2^(J+1)-1,:);
pvalue = sortrows(pvalue,1);

i = 1;
while (i <= m && pvalue(i,1) <= qm*i/m)
    i = i+1;
end

if i > 1
    for p = 1:i-1
        j = pvalue(p,2);
        k = pvalue(p,3);
        counts2mat = cell2mat(counts(:,j-j0+1));
        beta_jk = (2^(j/2))*(mean(counts2mat(:,2*k+1),1)-mean(counts2mat(:,2*k+2),1))/sqrt(T);
        rate_est = rate_est + 2^(j/2)*haar_func((2^j)*t-k*T,1)*beta_jk;
    end
end


end

