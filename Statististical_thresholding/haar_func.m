function [ psi ] = haar_func( t,T )

psi = ((t >= 0 & t<T/2) - (t >= T/2 & t<T))/sqrt(T);


end

