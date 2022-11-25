function g = bin2gray(b)
g(1) = b(1);
for i = 2 : length(b);
    x = xor(b(i-1), b(i));
    g(i) = x;
end
