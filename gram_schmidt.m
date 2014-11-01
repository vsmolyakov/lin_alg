%A=QR
%Q'Q=I and R is upper triangular
%othogonal-triangular decomposition

clc;
%m=8; n=8;
%randn('seed',0);
%A=randn(m,n);

%numerical loss of orthogonality for low precision:
%A=[0.70000 0.70711; 0.70001 0.70711]; 

%% Numerical Stability
randn('seed',0);
[U,X]=qr(randn(80)); %set U to a random orthogonal matrix
randn('seed',1);
[V,X]=qr(randn(80)); %set V to a random orthogonal matrix

S=diag(2.^(-1:-1:-80)); %set S to badly conditioned signular values
A=U*S*V;                %construct A=USV'=QR=QDU
                        %expect: diag(S)=diag(D)=diag(sigma_1,...,sigma_n)
[m,n]=size(A);

QC=zeros(m,n);
RC=zeros(m,n);
%% Classical Gram-Schmidt (GS): numerically unstable, complexity O(n^2)
for j=1:n
    v=A(:,j);
    for i=1:j-1
        RC(i,j)=QC(:,i)'*A(:,j); %compute r_ij = <aj,qi>
        v=v-RC(i,j)*QC(:,i);     %compute v=q_j orthogonal to {q_1,...,q_j-1}
    end
    RC(j,j)=norm(v,2); %compute normalization constant for v=q_j
    QC(:,j)=v/RC(j,j); %compute q_j
end

%% Modified Gram-Schmidt (MGS), complexity O(n^2)
QM=zeros(m,n);
RM=zeros(m,n);
VM=A;
for i=1:n
    RM(i,i)=norm(VM(:,i),2);
    QM(:,i)=VM(:,i)/RM(i,i);
    for j=i+1:n
        RM(i,j)=QM(:,i)'*VM(:,j);
        VM(:,j)=VM(:,j)-RM(i,j)*QM(:,i);
    end
end

%RM
%QM

%% Matlab Gram-Schmidt: qr(A), Householder triangularization, complexity
[Q_m,R_m]=qr(A);

%RC
%QC

%% Comparison Metrics

%frobenuous norm
norm(RC-R_m,'fro')
norm(QC-Q_m,'fro')
norm(RM-R_m,'fro')
norm(QM-Q_m,'fro')

%singular value sensitivity: diagonal(R)=diag(S)
figure;
semilogy(diag(RC),'-r','linewidth',1.5); hold on;
semilogy(diag(RM),'-b','linewidth',1.5); hold on;
semilogy(2.^-(1:length(diag(RC))),'--m');
legend('Classical GS', 'Modified GS', 'Singular Values: 1/2^j');   grid on;
title('Plot of the Diagonal Entries of R');

%orthogonality sensitivity: QQ'=I
norm(QC*QC'-eye(size(QC)),'fro')
norm(QM*QM'-eye(size(QC)),'fro')


%% Hadamard Inequality: |det(A)|<= Prod ||Aj||_2 for j=1...n
abs_det=abs(det(A))
abs_det_ub=1;
for j=1:n
    abs_det_ub = abs_det_ub * norm(A(:,j),2);
end
abs_det_ub



