function  U = extract_bearings(Z,R)
    % compute the (3 x nE) baseline (or bearing) matrix
    
    % Z: matrixe X (4n x 4n) con trasf. pairwise
    % R: cell array di matrici R (3 x 3)  con le rotazioni assolute [optional]
    %
    % U is the 3 x nE containing the baseline directions as columns,
    % in the same order as the edges in find_edge(A)
    
    % la matrice A (n x n) di adiacenza del grafo
    % la ottengo "strizzando" la X, sommando i suoi blocchi 4x4
    A = spones(blocksum(Z,4));
   
    [I,J]=find_edge(A);
    nedges=length(I); 
    U = zeros(3, nedges);
    for k=1:nedges
        i=I(k); j=J(k);
        if nargin == 2 % if rotation provided
            U(:,k) = -R{i}'*Z(4*i-3:4*i-1,4*j);
        else
            U(:,k) = Z(4*i-3:4*i-1,4*j);
        end
    end
end
