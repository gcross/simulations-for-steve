function [A,s,B] = HPhi
% Get a factorized Hermitian form for HPhi

A = zeros(3,3,4);
B = zeros(5,5,4);
s = zeros(4,1);

k=1;
s(k)=3/4;
A(:,:,k) = diag([1 1 0]);
B(:,:,k) = diag([1 1 0 0 0]);

k=2;
s(k)=1/4;
A(:,:,k) = diag([1 -1 0]);
B(:,:,k) = diag([1 -1 0 0 0]);

k=3;
s(k)=1/4;
A(:,:,k) = accumarray({2 1},1,[3 3]) + accumarray({1 2},1,[3 3]);
B(:,:,k) = accumarray({1 2},1,[5 5]) + accumarray({2 1},1,[5 5]);

k=4;
s(k)=-1/4;
A(:,:,k) = i*( accumarray({1 2},1,[3 3]) - accumarray({2 1},1,[3 3]) );
B(:,:,k) = i*( accumarray({2 1},1,[5 5]) - accumarray({1 2},1,[5 5]) );

end