%Non-negative Matrix Factorization (NMF) for topic extraction
%min norm(A-WH,'fro'), s.t. W,H>=0 and rank(W)=rank(H)=K (num of topics)

%% synthetic dataset
rng('default');
[A] = generate_data();

%% NMF
Kmax=min(size(A)); num_topics=linspace(1,Kmax,Kmax);
Kmax=length(num_topics); RMSE=zeros(Kmax,1); nnmf_time=zeros(Kmax,1);
for k=1:Kmax
    [~,RMSE(k),nnmf_time(k)]=nnmf_topic(A,num_topics(k));
end

figure;
[ax]=plotyy(1:Kmax,RMSE,1:Kmax,nnmf_time);
ylabel(ax(1),'RMSE'); ylabel(ax(2),'time, s');
xlabel(ax(2),'number of clusters');
title('Topic Extraction Using NNMF on Synthetic Data');
legend('RMSE','time');

K=2; %for 2D visualization
[nnmf_params,~,~]=nnmf_topic(A,K); 

%normalize beta
nnmf_beta_norm=nnmf_params.beta./repmat(sum(nnmf_params.beta),size(nnmf_params.beta,1),1);

%% generate plots
%plot word distributions for each of the K topics
figure;
for topic=1:K
    subplot(1,K,topic);
    plot(nnmf_beta_norm(:,topic)); hold on; axis tight;
    title('Non-Negative Matrix Factorization'); xlabel('words'); ylabel('\beta_{w|z}')
    legend_str=['topic ',num2str(topic)]; legend(legend_str,'Location','northeast');
end

figure;
biplot(nnmf_params.beta,'scores',nnmf_params.theta');
xlabel('Topic 1'); ylabel('Topic 2'); %axis([0 1 0 1]);

