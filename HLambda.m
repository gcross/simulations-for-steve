function [A,s,B] = HLambda(L)
% Find a factorized Hermitian operator basis for HLambda

if L==0
    s=1;
    A(:,:,1) = diag([0 0 0 0 1]);
    B(:,:,1) = diag([0 0 1]);
    return
end

A = zeros(5,5,9);
B = zeros(3,3,9);
s = zeros(9,1);

k=1;
s(k) = 1;
A(:,:,k) = diag([0 0 0 0 1]);
B(:,:,k) = diag([0 0 1]);

k=2;
s(k) = L^2/2;
A(:,:,k) = diag([0 0 0 1 0]);
B(:,:,k) = diag([1 0 0]);

k=3;
s(k) = L^2/2;
A(:,:,k) = diag([0 0 1 0 0]);
B(:,:,k) = diag([0 1 0]);

k=4;
s(k) = -(L/sqrt(2))/2;
A(:,:,k) = accumarray({3 5},1,[5,5]) + accumarray({5 3},1,[5,5]);
B(:,:,k) = accumarray({2 3},1,[3,3]) + accumarray({3 2},1,[3,3]);

k=5;
s(k) = (L/sqrt(2))/2;
A(:,:,k) = i*( accumarray({3 5},1,[5,5]) - accumarray({5 3},1,[5,5]) );
B(:,:,k) = i*( accumarray({2 3},1,[3,3]) - accumarray({3 2},1,[3,3]) );

k=6;
s(k) = (L/sqrt(2))/2;
A(:,:,k) = accumarray({4 5},1,[5,5]) + accumarray({5 4},1,[5,5]);
B(:,:,k) = accumarray({3 1},1,[3,3]) + accumarray({1 3},1,[3,3]);

k=7;
s(k) = (L/sqrt(2))/2;
A(:,:,k) = i*( accumarray({4 5},1,[5,5]) - accumarray({5 4},1,[5,5]) );
B(:,:,k) = i*( accumarray({3 1},1,[3,3]) - accumarray({1 3},1,[3,3]) );


k=8;
s(k) = -(L^2/2)/2;
A(:,:,k) = accumarray({3 4},1,[5,5]) + accumarray({4 3},1,[5,5]);
B(:,:,k) = accumarray({2 1},1,[3,3]) + accumarray({1 2},1,[3,3]);

k=9;
s(k) = (L^2/2)/2;
A(:,:,k) = i*( accumarray({3 4},1,[5,5]) - accumarray({4 3},1,[5,5]) );
B(:,:,k) = i*( accumarray({2 1},1,[3,3]) - accumarray({1 2},1,[3,3]) );


s = s/(1+L^2);

end