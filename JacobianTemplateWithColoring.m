function [z, JTW, JPIs] = JacobianTemplateWithColoring(functions, X, W, JPIs, dependency, Extra)
% JacobianTemplate applies multi-level reverse mode to get Jacobian
% efficiently
%
% functions - column cell array of function names
% X - column cell array of the initial inputs
% dependency - N-by-2 cell matrix, first column specifies dependency on X
%                               second column specifies dependency on Y
% Extra - extra parameters to define the computation
% z - result of the computation
% J - Jacobian matrix
N_funs = size(functions,1);
assert( 1 == size(functions,2), 'functions must be a column cell of strings');
assert( N_funs == size(dependency,1), 'dependency of every function must be specified');

[N_x, nrow] = size(X);
assert( 1 == nrow, 'input X must be a column cell array');
[xDim, nrow] = size(cell2mat(X));
assert( 1 == nrow, 'each element of X must be a column vector');

Y = cell(N_funs,1);

for n = 1:N_funs
    input = [cell2mat(X(dependency{n,1})); cell2mat(Y(dependency{n,2}))];
    Y{n} = feval(functions{n}, input, Extra);
end
z = Y{end};
assert( size(W,1) == size(z,1), 'output dimension must equal that of the dual vector');

idx_X = cell(length(X),1);  tail = 0;
for n = 1:length(X)
    idx_X{n} = (tail+1):(tail+length(X{n}));
    tail = tail + length(X{n});
end

JTW = zeros(xDim,size(W,2));
G = cell(N_funs,1);
G{N_funs} = W;
if isempty(JPIs)
    JPIs = cell(N_funs,1);
end
for n = N_funs:-1:1
    input = [cell2mat(X(dependency{n,1})); cell2mat(Y(dependency{n,2}))];
    dim_input = length(input); dim_output = length(Y{n});
    
    if dim_output > 4
        if isempty(JPIs{n})
            [JPIs{n}, ~] = getjpi(functions{n}, dim_output, dim_input, Extra);
        end
        [~, Jn] = evalj(functions{n}, input, Extra, dim_output);
        w = Jn'*G{n};
    else
        AD_fun = ADfun(functions{n},size(Y{n},1));
        options = setopt('revprod',G{n});
        [~, w] = feval(AD_fun, input, Extra, options);
        if 1 == size(w,1)
            w = w';
        end
    end
    
    idx = reshape(cell2mat(idx_X(dependency{n,1}))',1,[]);
    nx = length(idx);
    if 0 ~= nx
        JTW(idx,:) = JTW(idx,:) + w(1:nx,:);
    end

    head = nx + 1;
    for k = dependency{n,2}
        tail = head + size(Y{k},1) - 1;
        if isempty(G{k})
            G{k} = w(head:tail,:);
        else
            G{k} = G{k} + w(head:tail,:);
        end
        head = tail + 1;
    end
end