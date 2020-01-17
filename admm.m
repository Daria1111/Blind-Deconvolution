clear; clc;
pic = checkerboard(1);
m = pic(:);

h = randn(64,1);
X = h * m';
x = diag(X);

A = randn(64,64);
b = A * x;
ans = tracelasso(A,b);
figure
imshow(pic)

res = ans./h;
res = reshape(res,[8,8]);
figure
imshow(res)
%res

norm(A * (x - ans), 'fro') / norm(b, 'fro')
