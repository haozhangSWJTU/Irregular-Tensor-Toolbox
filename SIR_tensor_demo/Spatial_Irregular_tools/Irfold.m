function  [X] =Irfold(XTi,Xtrue)

X = NaN(size(Xtrue));
Xslice = Xtrue(:,:,1);
oo2 = find(~isnan(Xslice));
for i = 1:size(Xtrue,3)
    temp = X(:,:,i);
    temp(oo2) = XTi(:,i);
    X(:,:,i) =  temp;
end

end