function D = euclideanDistance(X1,Y1,X2,Y2)
    % D = euclideanDistance(X1,Y1,X2,Y2)
    % X1,Y1,X2,Y2 must all be of the same dimensions.
    
    DX = X1-X2;
    DY = Y1-Y2;
    D = sqrt(DX.^2+DY.^2);
    
    
end