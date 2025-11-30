function  [XTi,oo2] =Irunfold(Xtrue)
    
Xslice = Xtrue(:,:,1);
oo2 = find(~isnan(Xslice));
for i = 1:size(Xtrue,3)
    temp = Xtrue(:,:,i);
    XTi(:,i) =   temp(oo2);
end

end