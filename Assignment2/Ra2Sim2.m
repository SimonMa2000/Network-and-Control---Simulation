% Reading Assignment 2 EE130 By Simon
% Implementation of algorithms in paper:
% Revisiting Finite-Time Distributed Algorithms via 
% Successive Nulling of Eigenvalues
% Sam Safavi, Student Member, IEEE, andUsmanA.Khan, Senior Member, IEEE
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6876198
% Example 2 in the paper
function Ra2sim2()
clear;clc;
%% Part 1
% give prameters
% W = mat[nNode, nNode]
firstrow = [0 1 0 0 0 1];
W=circulant(firstrow,1);
% set initial state of the system 
% x0 = [nNode,1]
x0=[1 2 3 4 5 6]';
%%
nNode = size(W,1);
% lambda=[nNode, 1]
[V, D] = eig(W)
% W*V-V*D=eps[nNode,nNode]
lambda = diag(D,0)
%%
% find the 1/sqrt(nNode)*ones(nNode) column vector in V
iones=findOnesVector(V);
% delete repeating terms lambda
eignum=findDistinctEig(lambda);
delItem = @(Z,A)(Z(find(Z~=A)) );
%%
% set a_i, a_j 
% lambdad - distinct lambda (eigenvalues)
lambdad=lambda(eignum); %[nDist,1]
kbig=length(lambdad);
% get the specific lambda that yields eigenvector parallel to ones
templambdad=lambdad;
% coefa=lambdad;
%% calculate coef a for different eigenvector
% find the all-one=eigenvector
klambdad=find(lambdad==lambda(iones));
% MODIFICATION
% find -1 so yields error
% but need try catch statement maybe...

if klambdad==[]
    error('no all one eigenvector');
end
% prepare to calculate coefa
for k=1:kbig
% templambdad(klambdad) specific lambda_k
templambdad=delItem(lambdad,lambdad(klambdad));
tempprod=prod(templambdad-lambdad(klambdad)*ones(kbig-1,1));% multiply all (aj-ai)
% coef a
coefa=lambdad;
coefa(k)=1/tempprod+lambdad(klambdad);
end
% coefa [kbig, kbig] coef calculated from lambdak
%% calculate  eq. (5) in the letter
% coef=[kbig,1]
AtMultiply=eye(nNode, nNode);
At=zeros(kbig, nNode, nNode);
for k=1:kbig
    attemp=coefa(k)*eye(nNode, nNode);
    At(k,:,:) = attemp-D %diag( coefa(k)*ones(1,nNode) )-D;
    temp=squeeze(At(k,:,:));
%coefa(:,4)
    AtMultiply=AtMultiply*squeeze(At(k,:,:))
%disp({'a',k});
%squeeze(At(k,:,:))
end
% show
disp('Rank of V*AtMultiply*V^T');
rankFinal=rank(V*AtMultiply*V')
xfinal=V*AtMultiply*V'*x0;
disp('Final state achieved');
% lbnew=[-2 0 0 0 2]'
% lambda correspondance v1->2, v2->-2, v3->0
eigenvectorSpace= V(:,eignum)
% For this case, eigenvectorSpace= V(:,[1 2 4 6])
% This is bound to be since it yields to all one vector.
% as long as col 6 is in the subspace,it would be zero.
% angle between vector x4 and subspace of all eigenvectors 
angle=subspace(xfinal,eigenvectorSpace);
if (angle<1e-15)
disp('As we see, the result of iteration is in the subspace that consist of all eigenvectors')
end
% coefa
end
%%
function iones=findOnesVector(V)
 
% function iones=findOnesVector(V)
% V =[nNode, nNode] eigenvector of W
% return iones [1*1] column number of the first ones eigenvector  
col=size(V,2);
nNode=col;
iones=-1;
    for k=1:col
        if ( sum( abs(V(:,k)-1/sqrt(nNode)*ones(nNode)))<1e-8)
            disp('Good! found ones vector!'); 
            iones=k;
            break;
        end
    end
    if(iones==-1)
        disp('No all-one-vector!')
    end
end

%%
function eignum=findDistinctEig(lambda)

% eignum=findDistinctEigvalue(lambda)
% input lambda=diag(D,0),where D is  sorted eigenvalues of the given matrix
% output eignum=[1,nSubnum] subnumber of each distinct eigenvalue
% k flag for eignum k
% m flag for lambda [bigk,1]
% set first 
    k=1; m=1;    
    lambdanorep(m)=lambda(k);
    eignum(1)=1;
    % prep while loop
    bigk=length(lambda); 
    while k<bigk
        k=k+1;
        % if not equal add to eignum
        if (lambdanorep(m)~=lambda(k))
            m=m+1;
            lambdanorep(m)=lambda(k);
             % Track number of k which add in loop
            eignum(m)=k;
        end
    end
end
