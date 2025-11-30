function Z_filled = interpolate_nans(Z)

    [rows, cols] = size(Z);
    [X, Y] = meshgrid(1:cols, 1:rows);
    validIdx = ~isnan(Z);
    X_valid = X(validIdx);
    Y_valid = Y(validIdx);
    Z_valid = Z(validIdx);
    Z_filled = griddata(X_valid, Y_valid, Z_valid, X, Y, 'linear');
    nanIdx = isnan(Z_filled);
    if any(nanIdx(:))
        Z_filled(nanIdx) = griddata(X_valid, Y_valid, Z_valid, X(nanIdx), Y(nanIdx), 'nearest');
    end   
    Z_filled = reshape(Z_filled,[size(Z,1),size(Z,2),size(Z,3)]);
end