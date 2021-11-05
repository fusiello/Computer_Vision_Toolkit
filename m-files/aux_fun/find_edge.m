
function [out1, out2] = find_edge (A,in1,in2)
    % FIND_EDGE Locate edge in graph represented by A
    %
    % A is the nv x nv (symmetric) adjacency matrix associated to a
    % graph:
    % A(i,j)=1 if the edge (i,j) appears in the graph,
    % A(i,j)=0 otherwise.
    %
    % Syntax (same as MATLAB findedge):
    %
    % [sOut,tOut] = find_edge(A)
    % [sOut,tOut] = find_edge(A,idx)
    % idxOut = find_edge(A,s,t)
    %
    % Description
    %
    % [sOut,tOut] = find_edge(A) returns the source and target node
    % IDs, sOut and tOut, for all of the edges in graph G.
    %
    % [sOut,tOut] = find_edge(A,idx) finds the source and target nodes
    % of the edges specified by idx.
    %
    % idxOut = find_edge(A,s,t) returns the numeric edge index,
    % idxOut, for the SINGLE edge specified by the source and target node
    % pair s and t. Edge indices ranges from 1 to nE, where nE is the
    % number of edges in the graph. An edge index of 0 indicates an
    % edge that is not in the graph. The edges are ordered using
    % Matlab's convention on sparse matrices.
    %
    %
    % Authors: Federica Arrigoni, Beatrice Rossi, Andrea Fusiello, 2015
    
    % Consider the lower-triangular adjacency matrix, i.e. the edges
    % (i,j) and (j,i) are considered as a single edge.
    
    
    [I,J]=find(tril(A,-1));
    
    if nargin == 3  && nargout == 1 %  idxOut = find_edge(A,s,t)
        
        if length(in1) > 1 || length(in2) > 1
            error('s and t must be scalars')
        end
        
       s = in1; t = in2;
        
        % Invert t and s if t>s
        if t>s
            temp=t; t=s; s=temp;
        end
        
        % Find the index corresponding to the given edge
        % [~,out1]=ismember([s(:),t(:)],[I,J],'rows');
        out1 = find(I == s(:) & J == t(:));
        % out1 =  idxout;
        
    elseif nargin < 3 && nargout == 2
        %  [sOut,tOut] = find_edge(A)
        %  [sOut,tOut] = find_edge(A,idx)
        
        if nargin == 2
            idx = in1;
            out1 = I(idx); % sOut
            out2 = J(idx); % tOut
        else
            out1 = I; % sOut
            out2 = J; % tOut
        end
        
    else
        error('Invalid format');
    end
    
end





