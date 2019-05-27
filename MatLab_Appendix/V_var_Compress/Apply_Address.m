function ans = Apply_Address(ad,SIZE)

n = length(ad);

A = [[0.5,0];[0,0.5]];

x = [0,0]' ; 

b1 = [0,0]' ; 
b2 = [0,0.5]' ; 
b3 = [0.5,0]' ; 
b4 = [0.5,0.5]' ;


for i = 1:n 
    j = ad(n+1-i);
    
    if j == 1
        x = A*x + b1;
        
    elseif j ==2
        x = A*x + b2;
        
    elseif j ==3
        x = A*x + b3;
        
    elseif j ==4
        x = A*x + b4;
        
        
    end
end    
    ans = x*SIZE + 1; 
        