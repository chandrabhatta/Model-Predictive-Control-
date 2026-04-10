function z = lunarTerrain(x, y)

    % sloped plane
    ax = 0.05;  % slope along x
    ay = 0.02;  % slope along y
    c  = 0.1;   % offset
    z = ax * x + ay * y + c;

    % shallow craters
    crater_centers = [3.6 0.24; 0.37 1.84];
    crater_depths = [0.08, 0.06];
    crater_sizes  = [0.4, 0.3];

    for i = 1:length(crater_depths)
        xc = crater_centers(i,1);
        yc = crater_centers(i,2);
        A  = crater_depths(i);
        s  = crater_sizes(i);

        z = z - A * exp(-((x-xc).^2 + (y-yc).^2)/(2*s^2));
    end

end