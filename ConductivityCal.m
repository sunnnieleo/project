function [Map] = ConductivityCal( positions, particle, numPSystem)
%ConductivityCal Calculated conductivuty for a region y determining the
%number of particle reminaing in the region
%   Region size fixed (75X50)
%   Type of material must be defined as a character n, p, or d for n-doped
%   p-doped and depletion region
%

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin
ymax = 50e-9;
mass = 9.109E-31; %in kg
TauMN = 0.2E-12;%mean skattering time
nSum =0;
dSum =0;
pSum =0;
xmax = 75;
ymax = 50;
n_d = 25;
d_p = 50;

%define thermal velocity (source:
%https://en.wikipedia.org/wiki/Thermal_velocity)

for i=1:numPSystem
    %determine if particle in this region
    if positions (i, 1) <=n_d
        nSum = nSum+1;
        
    elseif (positions (i, 1) > n_d) && (positions (i, 1) <= d_p)
        dSum =dSum+1;
        
    else %positions (i, 1) >=d_p
        pSum = pSum+1;
        
    end
end

if particle == 'electron'
    charge = -1602E-19;
elseif particle == 'hole'
    charge = 1602E-19;
else
    error('invalid particle type')
end

conductivtyN = charge^2*TauMN*nSum/mass;
conductivtyD = charge^2*TauMN*dSum/mass;
conductivtyP = charge^2*TauMN*pSum/mass;


for i = 1:xmax
    for j = 1:ymax
        if i <=n_d
            Map(i, j) = conductivtyN;            
        elseif i > n_d && i < d_p
            Map(i, j) = conductivtyD;
        else %positions (i, 1) >=d_p
            Map(i, j) = conductivtyP;            
        end
    end
end

end

