function ans = Decode(C , Q)

Ind = code(9);

img = zeros(512,512);

[m,n] = size(Ind);

for i = 1:m
    for j = 1:n
        
        ad = char(Ind(i,j));
        
        img(i,j) = Apply_Ad(ad,C,Q);
        
    end 
end
ans = uint8(img);