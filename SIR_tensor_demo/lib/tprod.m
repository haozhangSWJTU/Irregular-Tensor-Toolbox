function C = tprod(A, B)
[n1, ~, n3] = size(A);
l = size(B, 2);
Af = fft(A, [], 3);
Bf = fft(B, [], 3);
Cf = zeros(n1, l, n3);
for i = 1 : n3
    Cf(:, :, i) = Af(:, :, i) * Bf(:, :, i);
end
C = ifft(Cf, [], 3);
end
