function [D,K,T,cnum]=descend_all(Parent,i)

% Find indexes of descendants of vertexes i from a tree 
% defined by the parent pointer Parent
% 
% Usage
%     [D,K,T]=descend_all(Parent,i)
%
% Inputs:
%     Parent is the parent pointer tree with 0 for the root
%     i is the index of the subtree root
% Outputs:
%     D is the list of descendent elements
%     K is the generation index of descendants, starting from i
%     T is the subtree with root i
%

% See also pred_all.m

if length(i)>1
    error('i must be length 1')
end

D=[i]; % indexes of descendants
K=[0]; % descendant rank (0 for the element itself)
J1=[i]; % elements to be treated at the next step
P = [1]; % position of elementf from J1 within the tree
T = [0]; % subtree descendant from i
cnum = []; % number of children

k=1; % descendant rank
while ~isempty(J1)
    L=length(T);
    J1_new=[];
    for i=1:length(J1)
        I=find(Parent==J1(i));
        J1_new=[J1_new I];
        T = [T; P(i)*ones(length(I),1)]; % adds children of J1(i) to T
        cnum(P(i)) = length(I); 
    end
    J1=J1_new;
    P = L+[1:length(J1)];
    D=[D J1];
    K=[K ones(size(J1))*k]; 
    k=k+1;
end

