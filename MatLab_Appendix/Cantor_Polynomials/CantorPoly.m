n_poly = 10;
Cantor = OrthogonalPolynomials(n_poly,[1/3,1/3],[0,2/3],[0.5,0.5],20);

for i = 1:n_poly
    %y = roots(sym2poly(Cantor{i}))
    figure(i)
    set(figure(i),'PaperSize',[6 10])
    
    fig = fplot(Cantor{i},[0,1]);
    
    
    set(gca,'YTick', [])
    
 
    saveas(fig,string(i)+'.png');
    %for j = 1:i
    %    vline(y(j),'r',string(y(j)))
    %end 
    %hline(0,'k')

    %hold on
    %coeffs(Cantor{i})
    %hold on
end
%hold off 
%set(gca, 'YScale', 'log')