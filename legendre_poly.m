%legendre polynomials
%orthogonal basis polynomials
%see also: lagrangian interpolation polynomials
 

%% legendre polynomials (cols of Q)

x=linspace(-128,128,128)'; %interval [-128,128]
A=[x.^0 x.^1 x.^2 x.^3];   %vandermonde matrix
[Q,R]=qr(A,0);

%normalize s.t. P1(x)=1
scale=Q(end,:);
Q=Q*diag(1./scale);

%plot(Q); axis tight; grid on; %plot cols of Q
plot(x,Q(:,1),'-r','linewidth',1.5); hold on;
plot(x,Q(:,2),'-g','linewidth',1.5); hold on;
plot(x,Q(:,3),'-b','linewidth',1.5); hold on;
plot(x,Q(:,4),'-m','linewidth',1.5); hold on;
axis tight; grid on; 
title('legendre polynomials');
xlabel('interval'); ylabel('P_j(x)');






