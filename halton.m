% van der Corput generator for Halton sequences
% b - base prime number, n - integer number recorded in a base b
% h(n,b) = sum_{k=0}^L*[d_k* b^{-k-1}]
%tic; for i=1:10000 a(i) = (halton(i,7)); end;for i=1:10000 b(i) = (halton(i,13)); end;figure; plot(a,b,'o');toc
%tic;figure;plot(rand(10000,1),rand(10000,1),'o');toc;
%figure; plot(a,b,'o')
function h = halton(n,b)
n_old = n;
h = 0;
f = 1/b;
while (n_old > 0)
   n_new = floor(n_old/b);
   r = n_old - n_new*b;
   h = h+f*r;
   f = f/b;
   n_old = n_new;
end
