clear;clc;


pic = checkerboard(1);


y = pic(:);


N = 64; 

L = 1;
B = randn(N,L);

alpha = randn(L,1);              
h = B*alpha; 


y = fft(y);
y_obs = diag(h) * y;

lambda = 1;
z = inv(diag(h));


x = z * y_obs;
ans = rpca(x,diag(h),y_obs,lambda);
n = norm(y_obs - diag(h) * ans) / norm(y_obs)
ans = ifft(ans);

ans = reshape(ans, [8,8]);
figure
imshow(abs(ans))


