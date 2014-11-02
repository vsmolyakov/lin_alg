%% Stable Algorithm
% ||fl(x)-f(fl(x))||/||f(fl(x))|| = O(eps_machine)
% for some fl(x) s.t. ||fl(x)-x||/||x|| = O(eps_machine), i.e.
% A stable algorithm gives nearly the right answer to nearly the right 
% question

% Backward Stability:
% Same as above, except fl(f(x))=f(fl(x)), i.e.
% A backward stable algorithm gives exact answer to nearly the right
% question

% Theorem: For a backward stable algorithm: 
% ||fl(f(x))-f(x)||/||f(x)|| = O(k(x)eps_machine), i.e.
% upper bound on accuracy increase linearly with stability number k(x)
% high numerical stability for small k(x) and vice versa

% Can estimate k(x) using global backward error stability rather than
% local forward error stability computations

% Householder triangularization with column pivoting is stable for
% almost all problems


%% Polynomial Example

%p(x)=(x-2)^9
coeff=[1 -18 144 -672 2016 -4032 5376 -4608 +2304 -512]';
x=linspace(1.920, 2.080, 160)';

p1=polyval(coeff,x);
p2=(x-2).^9;

figure;
plot(x,p1,'-r'); hold on;
plot(x,p2,'-b'); axis tight; grid on;
legend('fixed coeffs','recomputed coeffs');
title('sensitivity of (x-2)^9 to perturbations in x');
xlabel('perturbation of x'); ylabel('(x-2)^9');

%% Perturbation of Identify

r=1e-5;
m=2; n=2;
A=eye(m,n);
[V1,D1]=eig(A);
A(1,1)=A(1,1)+r;
[V2,D2]=eig(A);

%% Householder Triangularization (backward stable, eps_machine=1.11e-16)
R=triu(randn(50));   %R is upper triangular
[Q,X]=qr(randn(50)); %Q is unitary
A=Q*R;
[Q2,R2]=qr(A); %compute Q2 and R2 by Householder triangularization
norm(Q2-Q)/norm(Q)  %lost 13 digits of accuracy
norm(R2-R)/norm(R)  %lost 13 digits of accuracy
norm(A-Q2*R2)/norm(A)  %full eps_machine accuracy

%errors in Q2*R2 are highly correlated!
%to check: consider perturbing Q and R:
Q3=Q+1e-4*randn(50); %set Q3 s.t. it's closer to Q than Q2
R3=R+1e-4*randn(50); %set R3 s.t. it's closer to R than R2
norm(A-Q3*R3)/norm(A) %lost 12 digits of accuracy


%% SVD decomposition

m=50; n=50;
%randn('seed',1);

%synthesize A
[U,R1]=qr(randn(m,n));
[V,R2]=qr(randn(m,n));
S=diag(sort(rand(m,1), 'descend')); % uniform \in [0,1]
A=U*S*V';

%compute SVD
[U2,S2,V2]=svd(A);

norm(U-U2)
norm(V-V2)
norm(S-S2)
norm(A-U2*S2*V2')

plot(diag(U2*U'), '-r'); hold on;
plot(diag(V2*V'), '-b'); grid on;

cond(A)

%% Least Squares (regression)

m=100; n=15; %over-determined
t=(0:m-1)'/(m-1);
A=[];
for i=1:n A=[A t.^(i-1)]; end  %Vandermonde matrix
b=exp(sin(4*t));
b=b/2006.787453080206;   %normalization constant to make a15=1 for a15x^15

%Solve LS numerically and compute perturbation parameters
x=A\b; y=A*x;
kappa=cond(A) %
theta=asin(norm(b-y)/norm(b)) %angle between b and range(A)
eta=norm(A)*norm(x)/norm(y)   %proportion of length of y

%% Householder Triangularization

[Q,R]=qr(A,0);
x=R\(Q'*b);
x(end)

%solve by augmentation
[Q,R]=qr([A b],0);
Qb=R(1:n,n+1) %extract Q^h b
R=R(1:n,1:n)  %extract R
X=R\Qb;
x(end)

%solve via built-in Householder operator
x=A\b;
x(end)

%% Matrix Inverse (with SVD)

m=8; n=4;
A=randn(m,n);
[U,S,V]=svd(A);
S=diag(S);
tol=max(size(A))*S(1)*eps; %eps: machine eps
r=sum(S>tol); %number of singular values above tolerance level
S=diag(ones(r,1)./S(1:r));
X=V(:,1:r)*S*U(:,1:r)'; %X=inv(A)



