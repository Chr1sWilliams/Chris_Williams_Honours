function [M] = MomentMatrix(Alpha,Beta,P,n,it)
M = eye(n);
A = eye(n);
m = length(Alpha);
A = cell(1:m);

for k = 1:m
    a = eye(n);
    for i = 1:n
        for j = 1:n
            if j >= i 
                a(i,j) = nchoosek(j-1,i-1)* Alpha(k)^(i-1)*Beta(k)^(j-i) ;
            end
        end 
    end
    A{k} = a;
end

for i = 1:it
    M_sum = zeros(n);
    for k = 1:m
        M_sum = P(k).*ctranspose(A{k})*M*A{k}+ M_sum;
    end
    M = M_sum;
end
