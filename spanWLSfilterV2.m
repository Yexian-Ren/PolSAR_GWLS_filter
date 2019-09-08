function filterdata = spanWLSfilterV2(T3mat,lambda,alpha)

IN = T3mat(:,:,1)+T3mat(:,:,6)+T3mat(:,:,9);

[r,c] = size(IN);
k = r*c;

L=IN;

span1 = IN./(2*mean2(IN));
L(find(span1<=1)) = log(span1(find(span1<=1)));
L(find(span1>1)) = (span1(find(span1>1)))-1;

smallNum = 0.0001;

% Compute affinities between adjacent pixels based on gradients of L
dy = diff(L, 1, 1);
% dy = padarray(dy, [1 0], 'post');
dy = -lambda(1:end-1,:)./(abs(dy).^alpha(1:end-1,:) + smallNum);
dy = padarray(dy, [1 0], 'post');
dy = dy(:);

%×óÓÒ
dx = diff(L, 1, 2); 
% dx = padarray(dx, [0 1], 'post');
dx = -lambda(:,1:end-1)./(abs(dx).^alpha(:,1:end-1) + smallNum);
dx = padarray(dx, [0 1], 'post');
dx = dx(:);


% Construct a five-point spatially inhomogeneous Laplacian matrix
B(:,1) = dx;
B(:,2) = dy;
d = [-r,-1];
A = spdiags(B,d,k,k);

e = dx;
w = padarray(dx, r, 'pre'); w = w(1:end-r);
s = dy;
n = padarray(dy, 1, 'pre'); n = n(1:end-1);

D = 1-(e+w+s+n);
A = A + A' + spdiags(D, 0, k, k);

% Solve
% x = waitbar(0,'processing');
% for ii = 1:9
%     waitbar(ii/9)
%     temp = T3mat(:,:,ii);
%     res_t = A\temp(:);
%     filterdata(:,:,ii) = reshape(res_t,r,c);   
% end
% close(x)

% Suggestions from Dr.Qin for improving the algorithm
temp = reshape(T3mat,k,9);
res_t = A\temp;
filterdata = reshape(res_t,r,c,9);

end