function flambda = trianglesine_fun(x,lambda0,T,V,nu,A)

lambdam = 1;
chi = 1/10;
a = (2^(V+1)*chi/T)*lambdam;
lambdab = lambdam*(2-chi)/2;
Its = T*lambdam;

jvect = 0:2^(V+1);
intervals = repmat(T,1,2^(V+1)+1).*jvect/2^(V+1);
[X,I]=meshgrid(x,intervals);
index = max(I.*(I <= X),[],1)*2^(V+1)/T;
m = mod(index,2);
s = (m==0)-(m==1);

flambda(:,1) = (lambdab+a*(s.*x-s.*(index+m)*T/2^(V+1))+lambdam*A*sin(2^(nu+1)*pi*x/T+1/T))/Its;
flambda = lambda0*flambda + lambda0;

end