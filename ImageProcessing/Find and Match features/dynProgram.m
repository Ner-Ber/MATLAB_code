function mask = dynProgram(strip1, strip2)
% dynProgram creates a mask in the dimensions of the inputed patches.
% patch1/2 - strips of the two image being stiched. must be in same
% dimensions
% mask - the mask for the certain strip inputed

    e = (strip1 - strip2).^2;
    n = size(strip1, 1);
    E = zeros(size(strip1));    % definig the cumulative matrix
    E(1,:) = e(1,:);

    % fill the cumulative error matrix
    for i = 2:n
        % create matrix of minimal values of 3 upper neighboors
        Emin = min([E(i-1,:) ; [E(i-1,2:end), inf] ; [inf, E(i-1,1:end-1)]]);
        E(i,:) = Emin + e(i,:);
    end

    % create the mask by finding the minimal cost path
    mask = zeros(size(strip1));
    mask = padarray(mask, [0 1 0], inf);
    E = padarray(E, [0 1 0], inf);
    for i = n:-1:1
        if i == n
            [~, I] = min(E(i,:));
            mask(i,2:I) = 1;
        else
            [~, cor] = min(E(i,I-1:I+1));
            cor = cor-2;
            I = I+cor;
            mask(i,2:I) = 1;
        end
    end

    mask = mask(:,2:(end-1));


end
