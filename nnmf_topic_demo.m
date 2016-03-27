function [] = nnmf_topic_demo()
%Non-negative Matrix Factorization (NMF) for topic extraction
%min norm(A-WH,'fro'), s.t. W,H>=0 and rank(W)=rank(H)=K (num of topics)

%% synthetic dataset
rng('default');
[A] = generate_data();

figure; spy(A);
title('term document matrix'); xlabel('documents'); ylabel('word counts')

%% NMF
K = 4; %fix K
[nnmf_params,~,~]=nnmf_topic(A,K);

%% generate plots
%plot word distributions for each of the K topics
figure;
legend_str = {};
for topic=1:K
    stem(nnmf_params.beta(:,topic)); hold on; axis tight;
    title('Topic Matrix'); xlabel('words'); ylabel('\beta_{w|z}')
    legend_str{topic}=['topic ',num2str(topic)];
end
legend(legend_str,'Location','northeast');

%plot topic proportions for document matching
figure; imagesc(nnmf_params.theta); colormap(gray)
title('Topic Proportions'); xlabel('documents'); ylabel('topics');

end

function [nnmf_params,RMSE,time] = nnmf_topic(A,K)

%set zero entries to eps
zero_idx=find(A==0); A(zero_idx)=eps;

tic;
[W,H,RMSE]=nnmf(A,K);
time=toc;

%normalize
nnmf_params.beta=bsxfun(@times, W, 1./sum(W,1)); 
nnmf_params.theta=bsxfun(@times, H, 1./sum(H,1));

end

function [data] = generate_data()

N = 1e2;        %words (dictionary size)
D = 50;         %documents
K = randi(5)+2; %topics

data = zeros(N,D);
for k=1:K
    pi = dirrnd(ones(D,1)*0.05);    
    
    istart = floor(N/K*(k-1)) + 1;
    istop = floor(N/K*k);    
    
    numWords = randi(N,istop-istart+1,1);
    data(istart:istop,:) = mnrnd(numWords,pi);
end
%data = data';
data = sparse(data);
end

function pi = dirrnd(aa);
% pi = dirrnd(aa)
% draws a sample from a dirichlet with parameter vector aa

pi = randg(aa);
pi = pi/sum(pi);

end
