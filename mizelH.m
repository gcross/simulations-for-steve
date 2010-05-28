function H = mizelH(n,L)
% Builds the Hamiltonian for the 1-d Mizel model

% define site dimensions dq, db.
% blocks are (5)x(3x5)x(3x5)x...x(3x5)x(3x2)
dq = 5; db = 3;

M = 16*n; % total number of factored terms in H
N = 2*n+1; % total number of subsystems

% initialize H
H = init(M,N);

% define H_in
Hin = diag( [0 1 0 0 0] );

% Define all the on-site unitary terms with an identity gate
HU = kron( [1 -1;-1 1]/2 , eye(2) );  HU(dq,dq) = 0;

% add the initial term
H{1,1} = Hin+HU;

m=1; % initialize term counting (after the first term was added)

% H_U on-site terms (U = identity)
for j=3:2:(N-2)
    m=m+1;
    H{m,j} = HU;
end

% Bulk Phi term
[A,s,B] = HPhi;
r=length(s);
for j=2:2:(N-3)
    for k=1:r
        m=m+1;
        H{m,j} = s(k)*A(:,:,k);  H{m,j+1} = B(:,:,k);
    end
end

% Boundary Phi term
[A,s,B] = HPhin;
r=length(s);
j = N-1;
for k=1:r
	m=m+1;
	H{m,j} = s(k)*A(:,:,k);  H{m,j+1} = B(:,:,k);
end

% H_Lambda term
[A,s,B] = HLambda(L);
r=length(s);
for j=1:2:(N-2)
    for k=1:r
        m=m+1;
        H{m,j} = s(k)*A(:,:,k);  H{m,j+1} = B(:,:,k);
    end
end

% H_SF terms
for j=1:2:(N-2)
	m=m+1;
	H{m,j} = diag([1 1 1 1 0]);  H{m,j+1} = diag([0 0 1]);
	m=m+1;
	H{m,j} = diag([0 0 0 0 1]);  H{m,j+1} = diag([1 1 0]);
end


H = H(1:m,:); % delete any extra identity terms. (use "cell indexing")

end

function H = init(M,N)
% initialize the cell array for storing H

dq=5; db=3;

H = cell(M,N);
for m=1:M
	for j=1:N
        if j==N
            H{m,j}=eye(2);
        elseif mod(j,2)
            H{m,j}=eye(dq);
        else
            H{m,j}=eye(db);
        end
	end
end
% This is over allocating a bit, but we
% will componsate at the end by deleting the extra space.

end