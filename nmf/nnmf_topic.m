function [nnmf_params,RMSE,time] = nnmf_topic(A,K)

%set zero entries to eps
zero_idx=find(A==0); A(zero_idx)=eps;

%normalize A
A_norm=A./repmat(sum(A),size(A,1),1);

tic;
[W,H,RMSE]=nnmf(A_norm,K);
time=toc;

nnmf_params.beta=W;
nnmf_params.theta=H;

