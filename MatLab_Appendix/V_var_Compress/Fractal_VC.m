function [compress,Q] = Fractal_VC(img,k)

V = 4^k;

[compress,C] = kmeans(double(initialDivide(img,k+1)),V,'start','sample');



for i = k+2:9
     C = Next_Divide(C);
     [id,C] = kmeans(C,V,'start','sample');
     if i == k+1
         compress = id;
     else
         compress = [compress, id];
     end 
end
Q = int32(C(:,1));