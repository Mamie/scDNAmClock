s = csvread("~/scratch60/rep1_svd.csv", 1, 0);
A = csvread("~/scratch60/rep_1_data.csv", 1, 0);
M = size(A, 1);
N = size(A, 2);
sigma = 0.1;
lambda = sigma * 10;
svThreshold = 1E-6;

z           =   s(2:end);
s           =   [s(1) 1];
Is          =   1;
while( ~isempty(z) ),
idx         =   find(abs(z - s(Is, 1)) < svThreshold );
if( isempty(idx) )
  s           =   [s; [z(1) 1]];
z(1)        =   [];
Is          =   Is + 1;
end
z(idx)      =   [];
s(Is, 2)    =   s(Is, 2) + numel(idx);
end
clear z

%       warns the user about using SURE with a non-simple, not-full
%       rank matrix
if( any( s(:, 1) < svThreshold ) ),
fprintf('   +   [SURE_SVT] Warning: argument might be rank-deficient.\n')
end
if( any( s(:, 2) > 1 ) ),
fprintf('   +   [SURE_SVT] Warning: argument might have repeated singular values.\n')
end

%       find singular values above the threshold
idx_p   =   ( s(:, 1) > lambda );

x  	=   0;
if( any(idx_p)  ),
x   =   x + sum( 0.5*s(idx_p, 2).*(s(idx_p, 2) + 1) );
x   =   x + sum( (abs(M-N)*s(idx_p, 2) + 0.5*s(idx_p, 2).*(s(idx_p, 2) - 1)).*(max(0, s(idx_p, 1) - lambda)./s(idx_p, 1)) );
end

D  	=   zeros(size(s, 1));
for Ik = 1:size(s, 1),
D(:, Ik)    =   s(Ik, 2)*s(:, 2).*s(:, 1).*max(0, s(:, 1) - lambda)./(s(:, 1).^2 - s(Ik, 1).^2);
end
D( isnan(D) | isinf(D) | abs(D) > 1E6 )     =   0;

R 	=   x + 2*sum(D(:));
R           =   -M*N*sigma^2 + sum(min(lambda^2, s(:,1).^2)) + 2*sigma^2*R;