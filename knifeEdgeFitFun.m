function [w0fit,Zcenterfit] = knifeEdgeFitFun(Z,V)

% sort the rows from smallest to largest
[Z,sortIndex] = sortrows(Z);
V = V(sortIndex);

% normalized the voltage signal
minSub = (V - min(V));
T = minSub/max(minSub);

% determine whether increasing or decreasing
if T(1) > T(end)
    
    direction = -1;
    
else
    
    direction = +1;
    
end

% find the center of the Gaussian
ZcenterInd = find(T > 0.5);

if direction == -1

    Zcenter = Z(ZcenterInd(end));
    
    % estimate the beam waist
    ZleftInd = find(T > 0.9);
    Zleft = Z(ZleftInd(end));
    ZrightInd = find(T < 0.1);
    Zright = Z(ZrightInd(1));
    w0 = (Zright - Zleft)/2;

else
    
    Zcenter = Z(ZcenterInd(1));
    
    % estimate the beam waist
    ZleftInd = find(T < 0.1);
    Zleft = Z(ZleftInd(end));
    ZrightInd = find(T > 0.9);
    Zright = Z(ZrightInd(1));
    w0 = abs((Zright - Zleft)/2);
    
end

% fit the curve
modelFun = @(p,x)( p(1) + p(2)*p(1)*(erf((sqrt(2)/p(3))*(x-p(4))))); % fitting function
startingVals = [0.5,direction,w0,Zcenter]; % inititalize the parameters
coefEsts = nlinfit(Z,T,modelFun,startingVals);

% outputs
w0fit = coefEsts(3);
Zcenterfit = coefEsts(4);