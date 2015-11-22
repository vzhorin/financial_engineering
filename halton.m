% Van der Corput generators for Halton sequences
% b - prime number, n - interger number
% h(n,b) = sum_{k=0}^L(d_k b^{-k-1})
function h=halton(n,b)
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
