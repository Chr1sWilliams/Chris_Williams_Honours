function ans = Apply_Ad(str,C,Q)

s = length(str);

[m,n] = size(C);

V = length(Q);

k = s-n; 


g = 0;

for i = 1:k
    g = g + 4^(k-i+1)*(  double(string(str(i)))  -1); 
end


r = C(g +double(string(str(k+1)))  ,1);



    
 for j = 2:n

     In = 4*(r - 1) + double(string(str(k+j)));

     r = C(In,j); 
 end     
ans = round(Q(r));     
        
