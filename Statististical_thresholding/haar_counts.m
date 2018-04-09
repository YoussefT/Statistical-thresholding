function [ counts ] = haar_counts( points,T,J )

counts = zeros(1,2^J);

for i=0:2^J-1
     
    a = T*i/2^J;
    b = T*(i+1)/2^J;
    xab = sum(points < b)- sum(points < a);
    counts(i+1) = xab;
    
end
    
end