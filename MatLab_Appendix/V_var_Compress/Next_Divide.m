function ans = Next_Divide(C)

[V,n] = size(C);

L = sqrt(n);

ans = [];

for  j = 1:V
    S = C(j,:);
    S_v = vec2mat(S,L)';
    a = S_v(1:end/2, 1:end/2);
    b = S_v(1:end/2, end/2+1:end);
    c = S_v(end/2+1:end, 1:end/2);
    d = S_v(end/2+1:end, end/2+1:end);
    temp = [a(:)' ; b(:)' ; c(:)' ; d(:)'];
    ans = [ans ; temp];

end


% for j = 1:m
%     for i = 1:4
%         if i==1 && j == 1
%             ans = C(j,(i-1)*L + 1 : i*L );
%             
%         else 
%             ans = [ans; C(j,(i-1)*L + 1 : i*L )];
%             
%         end 
%         
%     end
% end 
%[idx,C] = kmeans(double(y),2)