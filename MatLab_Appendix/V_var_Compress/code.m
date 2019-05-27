function ans = code(n)

A = [["1","2"];["3","4"]];

for i =1:n-1
    
    A1 = "1"+A;
    A2 = "2"+A;
    A3 = "3"+A;
    A4 = "4"+A;

    A = [[A1,A2];[A3,A4]];
    
end    
ans=A;
        