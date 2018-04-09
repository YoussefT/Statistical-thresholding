function [ phi ] = haar_base_father( t,tau,T,j0)

j = j0;
k = 0:2^j0-1 ;

[KM,JM,TM,TAUM] = ndgrid(k,j,t,tau);
phi = 2.^(JM/2).*haar_func_father((2.^(JM)).*(TM+TAUM)-KM*T,T);

end

