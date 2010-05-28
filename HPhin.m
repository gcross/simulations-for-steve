function [A,s,B] = HPhin
% Factorized Hermitian form for the boundary Phi term

A = zeros(3,3,4);
B = zeros(2,2,4);
s = zeros(4,1);

k=1;
s(k)=3/4;
A(:,:,k) = diag([1 1 0]);
B(:,:,k) = diag([1 1]);

k=2;
s(k)=1/4;
A(:,:,k) = diag([1 -1 0]);
B(:,:,k) = diag([1 -1]);

k=3;
s(k)=1/4;
A(:,:,k) = accumarray({2 1},1,[3 3]) + accumarray({1 2},1,[3 3]);
B(:,:,k) = accumarray({1 2},1,[2 2]) + accumarray({2 1},1,[2 2]);

k=4;
s(k)=-1/4;
A(:,:,k) = i*( accumarray({1 2},1,[3 3]) - accumarray({2 1},1,[3 3]) );
B(:,:,k) = i*( accumarray({2 1},1,[2 2]) - accumarray({1 2},1,[2 2]) );


end