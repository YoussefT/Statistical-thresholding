function [ boolean ] = check_min_count( point_data,T,Jmax )

if not(iscell(point_data))
    point_cell{1} = point_data;
else
    point_cell = reshape(point_data, max(size(point_data)),1);
end

M = size(point_cell,1);
counts = zeros(M,2^(Jmax+1));
for m = 1:M
    counts(m,:) = haar_counts(point_cell{m},T,Jmax+1);
end

min_count = min(sum(counts,1));
boolean = min_count < 50 && M < 50;

end

