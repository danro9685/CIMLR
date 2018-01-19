function [K1, K2] = Estimate_Number_of_Clusters_SIMLR(alldata, NUMC)

for i = 1:length(alldata)
    
    if i == 1
        D_Kernels = multipleK(alldata{i});
        distX = mean(D_Kernels,3);
        W =  max(max(distX)) - distX;
        W = Network_Diffusion(W,max(ceil(size(alldata{1},1)/20),10));
    else
        D_Kernels = cat(3, D_Kernels,multipleK(alldata{i}));
        %D_Kernels = multipleK(alldata{i});
        distX = mean(D_Kernels,3);
        W0 =  max(max(distX)) - distX;
        W = W + Network_Diffusion(W0,max(ceil(size(alldata{1},1)/20),10));
    end
end
% distX = mean(D_Kernels,3);
% W =  max(max(distX)) - distX;
% W = Network_Diffusion(W,max(ceil(size(alldata{1},1)/max(NUMC)/2),10));

% [clusts,best_group_index,Quality,Vr] = cluster_rotate(W-diag(diag(W)),NUMC,1);Quality(isnan(Quality))=1;
% [clusts,best_group_index,Quality_plus,Vr] = cluster_rotate(W-diag(diag(W)),NUMC+1,1);Quality_plus(isnan(Quality_plus))=1;
% [clusts,best_group_index,Quality_minus,Vr] = cluster_rotate(W-diag(diag(W)),NUMC-1,1);Quality_minus(isnan(Quality_minus))=1;Quality_minus = [1, Quality_minus];

[Quality] = Estimate_Number_of_Clusters_given_graph(W, NUMC);
[Quality_plus] = Estimate_Number_of_Clusters_given_graph(W, NUMC+1);
[Quality_minus] = Estimate_Number_of_Clusters_given_graph(W, NUMC-1);

K1 = 2*(1 + Quality) - (2 + Quality_plus + Quality_minus);
K2 = K1.*(NUMC+1)./(NUMC);
subplot(1,2,1)
plot(NUMC,K1,'b-s','LineWidth',4);
title('Relative Quality')
subplot(1,2,2)
plot(NUMC,K2,'r-o','LineWidth',4);
title('Adjusted Quality')

end




function [quality] = Estimate_Number_of_Clusters_given_graph(W, NUMC)
%%%This function estimates the number of clusters given the two huristics
%%%given in the supplementary materials of our nature method paper
%W is the similarity graph
%NUMC is a vector which contains the possible choices of number of
%clusters.







if nargin < 2
    NUMC = 2:5;
end

if min(NUMC)==1
    warning('Note that we always assume there are more than one cluster.');
    NUMC(NUMC<1) = [];
end
W = (W + W')/2;

if ~isempty(NUMC)
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
    % compute unnormalized Laplacian
    L = D - W;
    degs(degs == 0) = eps;
    % calculate D^(-1/2)
    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
    % calculate normalized Laplacian
    L = D * L * D;
    % compute the eigenvectors corresponding to the k smallest
    % eigenvalues
    [U, eigenvalue] = eig(L);
    eigenvalue  = diag(eigenvalue);
    [a,b] = sort((eigenvalue),'ascend');
    eigenvalue = (eigenvalue(b));
    U = U(:,b);
    eigengap = abs(diff(eigenvalue));
    %eigengap = eigengap(NUMC-1);
    %     for i = 1:length(NUMC)
    %         eigengap(i) = (1+eigengap(i))/(1+mean(eigenvalue(1:NUMC(i))));
    %     end
    for ck = NUMC
        Cindex = find(NUMC==ck);
        if ck == 1
            quality(Cindex) = sum(sum(diag(1./(U(:,1)+eps))*U(:,1)));
        else
            UU = U(:,1:ck);
            UU = UU./repmat(sqrt(sum(UU.^2,2))+eps,1,size(UU,2));
            [EigenvectorsDiscrete,EigenVectors ]=discretisation(UU);
            EigenVectors = EigenVectors.^2;
            [temp1,temp] = sort(EigenVectors,2, 'descend');
            quality(Cindex) = (1-eigenvalue(ck+1))/(1-eigenvalue(ck))*sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
            %quality(Cindex) = sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
        end
    end
    % %     NUMC = NUMC+1;
    % %
    % %     eigengap_plus = abs(diff(eigenvalue));
    % %     eigengap_plus = eigengap_plus(NUMC-1);
    % %     %     for i = 1:length(NUMC)
    % %     %         eigengap_plus(i) = (1+eigengap_plus(i))/(1+mean(eigenvalue(1:NUMC(i))));
    % %     %     end
    % %     for ck = NUMC
    % %         Cindex = find(NUMC==ck);
    % %         UU = U(:,1:ck);
    % %         UU = UU./repmat(sqrt(sum(UU.^2,2))+eps,1,size(UU,2));
    % %         [EigenvectorsDiscrete,EigenVectors ]=discretisation(UU);
    % %         EigenVectors = EigenVectors.^2;
    % %         [temp1,temp] = sort(EigenVectors,2, 'descend');
    % %         %quality_plus(Cindex) = (1-eigenvalue(ck+1))/(1-eigenvalue(ck))*sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
    % %         quality_plus(Cindex) = sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
    % %     end
    % %
    % %
    % %     NUMC = NUMC-1;
    % %
    % %
    % %
    % %     NUMC = NUMC-1;
    % %
    % %
    % %     for ck = NUMC
    % %
    % %         Cindex = find(NUMC==ck);
    % %         if ck == 1
    % %             quality_minus(Cindex) = sum(sum(diag(1./(U(:,1)+eps))*U(:,1)));
    % %             %quality_minus(Cindex) = size(U,1);
    % %         else
    % %             UU = U(:,1:ck);
    % %             UU = UU./repmat(sqrt(sum(UU.^2,2))+eps,1,size(UU,2));
    % %             [EigenvectorsDiscrete,EigenVectors ]=discretisation(UU);
    % %             EigenVectors = EigenVectors.^2;
    % %             [temp1,temp] = sort(EigenVectors,2, 'descend');
    % %             %quality_plus(Cindex) = (1-eigenvalue(ck+1))/(1-eigenvalue(ck))*sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
    % %             quality_minus(Cindex) = sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
    % %         end
    % %     end
    % %
    % %     NUMC = NUMC+1;
    % %
    % %     adjusted_quality = quality;
    % %     quality = (2+2*quality) ./ (2 + quality_plus+quality_minus);
    % %     %adjusted_quality = quality.*(NUMC+1)./(NUMC);
    % %     %eigengap = (eps+eigengap_plus) - (eps+eigengap);
    % %     %[tt1, t1] = sort(eigengap,'descend');K1 = NUMC(t1(1));K12 = NUMC(t1(2));
    % %
    % %     %[tt2, t2] = sort(quality,'ascend'); K1 = NUMC(t2(1));K2 = NUMC(t2(2));
    % %     figure
    % %     plot(NUMC,quality,'b-s','LineWidth',4);
    % %     title('Separation Cost')
    % %     figure
    % %     plot(NUMC,quality.*(NUMC+1)./(NUMC),'r-o','LineWidth',4);
    % %     title('Adjusted Separation Cost')
    
end


end

function D_Kernels = multipleK(x)


N = size(x,1);
KK = 0;
sigma = [2:-0.25:1];
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 10:2:30;
t=1;
for l = 1:length(allk)
    if allk(l) < (size(x,1)-1)
        TT=mean(T(:,2:(allk(l)+1)),2)+eps;
        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
        Sig=Sig.*(Sig>eps)+eps;
        for j = 1:length(sigma)
            W=normpdf(Diff,0,sigma(j)*Sig);
            Kernels(:,:,KK+t) = (W + W')/2;
            t = t+1;
        end
    end
end

for i = 1:size(Kernels,3)
    K = Kernels(:,:,i);
    k = 1./sqrt(diag(K)+eps);
    G = K.*(k*k');
    %G = Network_Diffusion(G,10);
    
    %G = K;
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    %D_Kernels(:,:,i) = G;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
    
end

end




%
function W = Network_Diffusion(A, K)
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = (TransitionFields(P));
[U,D] = eig(P);
d = real((diag(D))+eps);
alpha = 0.8;
beta = 2;
d = (1-alpha)*d./(1-alpha*d.^beta);
D = diag(real(d));
W = U*D*U';
W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
%W=D*(W);
W = (W+W')/2;

end


function [EigenvectorsDiscrete,EigenVectors]=discretisation(EigenVectors)
%

[n,k]=size(EigenVectors);

vm = sqrt(sum(EigenVectors.*EigenVectors,2));
EigenVectors = EigenVectors./repmat(vm+eps,1,k);

R=zeros(k);
% R(:,1)=EigenVectors(1+round(rand(1)*(n-1)),:)';
R(:,1)=EigenVectors(round(n/2),:)';
%R(:,1)=EigenVectors(n,:)';
c=zeros(n,1);
for j=2:k
    c=c+abs(EigenVectors*R(:,j-1));
    [minimum,i]=min(c);
    R(:,j)=EigenVectors(i,:)';
end

lastObjectiveValue=0;
exitLoop=0;
nbIterationsDiscretisation = 0;
nbIterationsDiscretisationMax = 20;%voir
while exitLoop== 0
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;
    EigenvectorsDiscrete = discretisationEigenVectorData(EigenVectors*R);
    [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors+eps,0);
    NcutValue=2*(n-trace(S));
    
    if abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end
end

function Y = discretisationEigenVectorData(EigenVector)
% Y = discretisationEigenVectorData(EigenVector)
%
% discretizes previously rotated eigenvectors in discretisation
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

[n,k]=size(EigenVector);


[Maximum,J]=max(EigenVector');

Y=sparse(1:n,J',1,n,k);
% Y = J';
end
