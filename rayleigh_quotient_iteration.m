
%Rayleigh Quotient (least squares eigenvalue approximation)
%l_min < x'Ax / x'x < l_max
%min_a ||Ax-ax||_2

%Power Iteration
%v(k)=Av(k-1)=A^k v(0), v(k) converges to max eigenvector

%Iterative methods vs Direct (non-iterative) methods
%Iterative methods can lead to lower complexity


%% Rayleigh Quotient (eigenvalue estimate from an eigenvector estimate)
m=4; n=4;
A=randn(m,n);
A=1/2*(A+A'); 
[V,D]=eig(A);     %sorted (ascending) for symmetric A
x_eig=V(:,end);   %max eigen-vector
l_eig=D(end,end); %max eigen-value

%raleigh quotient
l=(x_eig'*A*x_eig)/(x_eig'*x_eig)

%perturb x_eig
x_eig2=x_eig+1e-1*rand(1);
l2=(x_eig2'*A*x_eig2)/(x_eig2'*x_eig2)

norm(l-l2)/norm(l)

%% Power Iteration (Von Mises iteration) (eigenvector estimate from an eigenvalue estimate)
num_iter=10;
w=zeros(m,1);

v=randn(m,1);
v=v/norm(v);
l=zeros(num_iter,1);
for k=1:num_iter
    w=A*v;         %apply A
    v=w/norm(w);   %normalize
    l(k)=v'*A*v;   %rayleigh quotient
end

figure;
plot(1:num_iter,l,'linewidth',2); grid on; axis tight;
title('Power Iteration of Largest Eigenvalue');
xlabel('number of iterations'); ylabel('max eigenvalue');

%% Inverse Iteration (speed up of eigenvector estimate from an eigenvalue estimate)
%  eig(A) <-> lambda
%  eig(A-mu*I) <-> lambda - mu
%  eig(A-mu*I)^-1 <-> 1/(lambda-mu)
%  mu: initial guess of eigenvalue
%  (A-mu*I)^-1 is large in the direction of eigenvector
%  corresponding to the closest eigenvalue of mu
%  convergence depends on the ratio of closest eig to 2nd closest to mu
%  use mu as initial guess to approximate any eigenvector of A

num_iter=10;
w=zeros(m,1);

[V,D]=eig(A);
mu=D(2,2)+1e-1*randn(1); %proximal to lambda_2

v=randn(m,1);
v=v/norm(v);
l=zeros(num_iter,1);
verr=zeros(num_iter,1);
for k=1:num_iter
    w=inv(A-mu*eye(size(A)))*v;   %apply (A-mu*I)^-1 
    v=w/norm(w);                  %normalize
    l(k)=v'*A*v;                  %rayleigh quotient
    verr(k)=norm(V(:,2)-v)/norm(V(:,2));
end

figure;
subplot(211)
plot(1:num_iter,l,'linewidth',2); grid on; axis tight;
title('Power Iteration of Closest Eigenvalue');
xlabel('number of iterations'); ylabel('max eigenvalue');
subplot(212)
plot(1:num_iter,verr,'linewidth',2); grid on; axis tight;
title('l2 Relative Error of Eigenvector Estimate');
xlabel('number of iterations'); ylabel('relative error');

%% Rayleigh Quotient Iteration (e.g. On-line PCA) (combination of rayleigh quotient and inverse iteration)

%iterative refinement of eigenvector and eigenvalue estimate through
%application of rayleigh quotient (eigenvalue) and inverse iteration
%(eigenvector)

%spectacular convergence, each iteration triples the number of digits
%of accuracy

num_iter=10;
w=zeros(m,1);
l=zeros(num_iter,1); 
verr=zeros(num_iter,1);

%init
v=randn(m,1); v=v/norm(v);  %unit norm randn eigenvector
l(1)=v'*A*v;                %corresponding init eigenvalue
for k=2:num_iter            %iterate until convergence
    w=inv(A-l(k-1)*eye(size(A)))*v;  %apply inv(A-mu*I)
    v=w/norm(w);                     %normalize
    l(k)=v'*A*v;                     %rayleigh quotient
end

figure;
plot(1:num_iter,l,'linewidth',2); grid on; axis tight;
title('Power Iteration of Largest Eigenvalue');
xlabel('number of iterations'); ylabel('max eigenvalue');

