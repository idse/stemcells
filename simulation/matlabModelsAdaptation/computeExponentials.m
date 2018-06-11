f = @(t,L,delta) exp(-delta*t)-exp(-L*t);

t = 0:1:100;

for ii = 1:20
    L = 0.1*ii;
    sol(ii,:) = f(t,L,0.01);
end
   

plot(sol');