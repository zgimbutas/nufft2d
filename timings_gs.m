
tic
ms = 50 
mt = 54 

n1 = ms^2/4
n2 = mt^2

%ms = 2
%mt = 1

%n1 = 2
%n2 = 2

%ms = 256*2
%mt = 256*2

%n1 = 256*2
%n2 = 256*2

k1 = [-n1/2: (n1-1)/2];
k2 = [-n2/2: (n2-1)/2];
nj = n1*n2; 

j = k2+n2/2+1 + (k1.'+n1/2)*n2;
j = j(:);

xx = pi * cos(-pi*k1/nj);
yy = pi * cos(-pi*k2/nj);
[xj, yj] = meshgrid(xx, yy);
cj = complex( sin(pi*j/n1), cos(pi*j/n2));

toc

format shorte

iflag = +1;

disp('=== direct ===')
tic
fk0 = dirft2d1(nj,xj,yj,cj,iflag,ms,mt);
toc

tic
cj0 = dirft2d2(nj,xj,yj,iflag,ms,mt,fk0);
toc


disp('=== type 1 ===')
eps = 1e-2
tic
  fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-4
tic
  fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-8
tic
  fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-12
tic
  fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)


disp('=== type 2 ===')
eps = 1e-2
tic
  cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk0);
toc

%%abs_error=norm(cj0-cj1,2)
rel_error=norm(cj0-cj1,2)/norm(cj0,2)

eps = 1e-4
tic
  cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk0);
toc

%%abs_error=norm(cj0-cj1,2)
rel_error=norm(cj0-cj1,2)/norm(cj0,2)

eps = 1e-8
tic
  cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk0);
toc

%%abs_error=norm(cj0-cj1,2)
rel_error=norm(cj0-cj1,2)/norm(cj0,2)

eps = 1e-12
tic
  cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk0);
toc

%%abs_error=norm(cj0-cj1,2)
rel_error=norm(cj0-cj1,2)/norm(cj0,2)

stop
disp('=== direct ===')

[f1, f2] = meshgrid([-ms/2:(ms-1)/2],[-mt/2:(mt-1)/2]);
tic
fk0 = dirft2d3(nj,xj,yj,cj,iflag,ms*mt,f1,f2);
toc

disp('=== type 3 ===')
eps = 1e-2
tic
  fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,ms*mt,f1,f2);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-4
tic
  fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,ms*mt,f1,f2);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-8
tic
  fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,ms*mt,f1,f2);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)

eps = 1e-12
tic
  fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,ms*mt,f1,f2);
toc

%%abs_error=norm(fk0-fk1,2)
rel_error=norm(fk0-fk1,2)/norm(fk0,2)
