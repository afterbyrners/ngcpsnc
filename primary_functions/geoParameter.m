function [lambda, S, AR, c_bar, x_mgc, y_mgc] = geoParameter(b,cr,ct,sweep) % various parameters that need to be found over and over
    lambda = ct/cr; % Taper Ratio
    S = 0.5*b*cr*(1+lambda); % Surface area
    AR = b^2/S; % Aspect Ratio
    c_bar = 2/3 * cr * (1+lambda+lambda^2)/(1+lambda); % self explanetory
    y_mgc = b/6 * (1+2*lambda)/(1+lambda); % same with this 
    x_mgc = y_mgc*tan(sweep); % and this as well
end