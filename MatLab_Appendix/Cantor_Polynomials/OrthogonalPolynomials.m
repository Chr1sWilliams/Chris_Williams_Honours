function [Poly] = OrthogonalPolynomials(n,Alpha,Beta,P,it)
syms x;
Poly = cell(1,n);
for i = 2:n+1
    A = MomentMatrix(Alpha,Beta,P,i,it);
    A = sym(A);
    for j = 1:i
        A(i,j) = x^(j-1);
    end
   Poly{i-1} = det(A);
end

