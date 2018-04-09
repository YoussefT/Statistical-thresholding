function [ psi ] = haar_base( t,tau,T,j0,J)

j = j0:J;
k = 0:2^J-1 ;

[KM,JM,TM,TAUM] = ndgrid(k,j,t,tau);
psi = 2.^(JM/2).*haar_func((2.^(JM)).*(TM+TAUM)-KM*T,T);

end

