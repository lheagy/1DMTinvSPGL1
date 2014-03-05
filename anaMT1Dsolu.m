function Z1D = anaMT1Dsolu(LayMat,Freq)
% anaMT1Dsolu.m
% Analytic solution for a horizontally layered earth.
%
% Output:
%   Z - 1D MT impedance. Length of Freq. Is in units [mV/km/nT]
%
% Inputs:
%   LayMat - Matrix with of 2 columns and N rows, where N is the number of
%       layers. Column 1 contains depth to the top of each layer and column
%       2 contains the conductivity value for each layer. A Halfspace is
%       assumed at the bottom.
%   Freq - frequencies to run for.
%
% GKR November 2010

%% Define constants
mu0 = 4*pi*1e-7;
eps0 = 8.85*1e-12;
untconv = 1; %1/(4*pi*1e-4); % Convert Z[V/A] to Z[mV/km/nT] (as in EDI)


Z1D = zeros(length(Freq),1);

% Calculate the layer thikness
d=abs(diff(LayMat(:,1)));

% Loop through  all frequencies
for nrFreq=1:length(Freq)
    % Calculate the values for the bottom layer
    Zh = zeros(length(LayMat),1);
    
    Zh(end) = (mu0*2*pi*Freq(nrFreq))/sqrt(mu0*eps0*(2*pi*Freq(nrFreq))^2 - 1i*mu0*LayMat(end,2)*2*pi*Freq(nrFreq));
    % Loop through all the layers
    for nrLay=(length(LayMat)-1):-1:1
        % Wave number
        k =  sqrt(mu0*eps0*(2*pi*Freq(nrFreq))^2 - 1i*mu0*LayMat(nrLay,2)*2*pi*Freq(nrFreq));
        Z = (mu0*2*pi*Freq(nrFreq))/k;  
       
        % Calculate the value
        Zh(nrLay) = Z*((Zh(nrLay+1)+Z*tanh(1i*k*d(nrLay)))/(Z+Zh(nrLay+1)*tanh(1i*k*d(nrLay))));
    end
    
%     Z1D(nrFreq) = untconv*Zh(1); % Took out 8.5.2012, should output in [V/A]
    Z1D(nrFreq) = Zh(1);
end
