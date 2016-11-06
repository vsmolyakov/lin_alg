% Page-Rank Algorithm

%% generate input
rng('default');
n=500; %number of links
i = randi(n,1,0.4*n);  %"to" 
j = randi(n,1,0.4*n);  %"from"
G = sparse(i,j,1,n,n); %connectivity matrix

%% page-rank
p = 0.85;        %prob that a random walk follows a link
delta = (1-p)/n; %prob of choosing a random page that doesn't follow a link
c = sum(G,1);    %out-degree of a page
k = find(c~=0);  %pages with outgoing links
D = sparse(k,k,1./c(k),n,n);
e = ones(n,1);
I = speye(n,n);

%power method
x = (I-p*G*D)\e; %Ax = x, where A = pGD is the state transition matrix
x = x/sum(x);    %page rank


%% generate plots

figure; spy(G,20); title('connectivity matrix');
figure; bar(x); xlim([0, n]); title('page rank');


