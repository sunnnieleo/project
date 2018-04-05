function [conductivty, pos, vel] = regionInit( particle, numPSystem)
%regionInit Initalize different types of materials to statify n, p or
%depletion region types
%   Type of material must be defined as a character n, p for n-doped
%   p-doped and depletion region
%
%   num particle in the system must be provided, and be divisible by 10
%   region size are fixed for this lab

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin
ymax = 50e-9;
xmaxNano = 75e-9;
mass = 9.109E-31; %in kg
TauMN = 0.2E-12;%mean skattering time
v_th = sqrt(k*T/mass);
% Conc = 0.5*numPSystem;
numP = numPSystem; 

%define thermal velocity (source:
%https://en.wikipedia.org/wiki/Thermal_velocity)

if isequal(particle ,'electron')    
    charge = -1602E-19; %in C
    %         numP = Conc;
    [pos, vel] = randInit( numP, 0, xmaxNano/2, ymax, v_th );
    conductivty = charge^2*TauMN*numP/mass;
    
    
elseif  isequal(particle ,'hole')  
    charge = 1602E-19; %in C
    
    %             numP = Conc;
    [pos, vel] = randInit( numP, xmaxNano/2, xmaxNano, ymax, v_th );
    conductivty = charge^2*TauMN*numP/mass;
else
    error('invalid particle type')
end

end

