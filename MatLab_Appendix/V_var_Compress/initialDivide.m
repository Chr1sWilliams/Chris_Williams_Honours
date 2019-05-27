function ans = initialDivide(img,k)

n = unique(size(img));

V = 4^k;

b_l = n/sqrt(V) - 1;

per = permn([1,2,3,4],k);

for i = 1 : length(per)
    
    add = Apply_Address(per(i,:),n);
  
    sub = img(add(1):add(1)+b_l, add(2):add(2)+b_l);
    
    sub = sub(:)' ; 

    if i ==1 
        MAT = sub;
    else
        MAT = [MAT;sub];
    end
end
ans = MAT;


