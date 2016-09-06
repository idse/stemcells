%%

clear pi;
L = 1000;
n = 5;
m = 3;
A = 10;
alphax = 0;
alphay = 1;
[X, Y] = meshgrid(1:L,1:L);
test = A*sin(n*2*pi*X/L + alphax).*sin(m*2*pi*Y/L + alphay);

%%
%plot(test(1,:))
F = fft2(test)*(2/L)^2;
absF = abs(F);
max(absF(:))
[i,j] = find(F == max(max(absF(1:100,1:100))),1,'first');
atan(imag(F(i,j))/real(F(i,j)))
[i,j]-1

%%
L = 1024;
F = abs(fft2(f)*(2/L)^2);
F(1,1) = 0;
[m,n] = find(F == max(F(:,1)),1,'first');
m = m-1;
n = n-1;

%%
n = 1;
m = 3;
A = 10;
Lx = size(S4Stitched,2);
Ly = size(S4Stitched,1);
[X, Y] = meshgrid(1:Lx,1:Ly);
f = A*sin(n*2*pi*X/Lx).*sin(m*2*pi*Y/Ly);
imshow(f,[])

%%
F = abs(fft2(f)*(2/Lx)*(2/Ly));
%F(1,1) = 0;
[A, idx] = max(F(:));
[m,n] = find(F == A,1,'first');
m = m-1;
n = n-1;
[n m]
A

%%
f = S4Stitched;
F = abs(fft2(f)*(2/Lx)*(2/Ly));
F(1,1) = 0;
[A, idx] = max(F(:));
[m,n] = find(F == A,1,'first');
m = m-1;
n = n-1;
[n m]
A
% phase