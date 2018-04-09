function [ rate_est ] = thresholding_innovation( points,t,T,j0,J,alpha )

%%% Estimate initial projection on V_j0
rate_est = thresholding_linear(points,t,T,j0);


M = size(points,1);
counts = cell(M,J-j0+1);
for j=j0:J
    for m = 1:M
        counts{m,j-j0+1} = haar_counts(points{m},T,j+1);
    end
end
pvalues = zeros(J+1,2);
Q = J-j0+1;

%%% Get the pvalue for each test

for j=j0:J
    counts2mat = cell2mat(counts(:,j-j0+1));
    R = pairconstrate_lrtstat(counts2mat,2^j);
    pvalues(j+1) = chi2cdf(R,2^j,'upper');
    pvalues(j+1,2) = j;
end

pvalues = pvalues(j0+1:J+1,:);

%%% Find minimal index for Holm-Bonferonni correction

pvalues = sortrows(pvalues,1);
im = 1;
while (im <= Q && pvalues(im,1) <= alpha/(Q+1-im))
    im = im+1;
end
   
%%% Keep resolutions where the test is rejected

if im > 1
    for i=1:im-1
        L = pvalues(i,2);
        counts2mat = cell2mat(counts(:,L-j0+1));
        betas = zeros(2^L,1);
        for k = 1:2^L
            betas(k) = (2^(L/2))*(mean(counts2mat(:,2*k-1),1)-mean(counts2mat(:,2*k),1))/sqrt(T);         
        end
        S = repmat(squeeze(betas),1,length(t));
        H = squeeze(haar_base(t,0,T,L,L));
        if L == 0
            H = H';
        end
        w = S.*H;
        w = sum(w,1);
        rate_est = rate_est + w;
    end
end

end

