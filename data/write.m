function write(filepath, outpath)
    problem = load(filepath);
    G = problem.Problem.A;

    % make the graph square if necessary
    [x,y] = size(G);
    if x ~= y
        [I,J,S] = find(G);
        n = max(x,y);
        G = sparse(I,J,S,n,n);
    end

    % find the largest connected subgraph and remove any vertices not in it
    [S,C] = graphconncomp(G, 'Weak', true, 'directed', true);
    [modee,freq] = mode(C);
    freq

    if graphisdag(G)
        disp('acyclic')
    else
        disp('cyclic')
    end

    [rows,cols,vals] = find(G);
    fileID = fopen(outpath, 'w');
    formatSpec = '%d %d\n';
    for i=1:length(rows)
        row = rows(i);
        col = cols(i);
        val = vals(i);

        % if both vertices are in the largest subgraph, then use
        if C(row)==modee && C(col)==modee
            fprintf(fileID, formatSpec, row, col);
        end
    end
end