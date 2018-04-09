function [ rate_est ] = thresholding_hard( points,t,T,j0,J,c )

%%% Estimate initial projection on V_j0
rate_est = thresholding_linear(points,t,T,J);


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
    for k = 0:2^j-1
        beta_var = (2^j)*(mean(counts2mat(:,2*k+1)+counts2mat(:,2*k+2),1))/T;
        beta_jk = (2^(j/2))*(mean(counts2mat(:,2*k+1),1)-mean(counts2mat(:,2*k+2),1))/sqrt(T);
        if abs(beta_jk) >= c*sqrt(beta_var)
            betas_rec(j-j0+1,k+1) = beta_jk;
        end
    end        
end

betas_rec = repmat(betas_rec',1,1,length(t));

psi_mat = haar_base(t,0,T,j0,J); % set up the mother wavelets family.

rate_est = rate_est + squeeze(sum(sum(betas_rec.*psi_mat,2),1))'; % combine all previous results in the final estimation.


end

