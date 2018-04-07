%% Project Sonya Stuhec- Leonard 100963181

%define electron parameters
mass = 9.109E-31; %in kg
chargeN = -1.602E-19; %in C
chargeP = 1.602E-19; %in C
TauMN = 0.2E-12;

k = physconst('Boltzmann'); %Use of constants in matlab
T = 300; % temperature in Kalvin

%define thermal velocity (source:
%https://en.wikipedia.org/wiki/Thermal_velocity)
v_th = sqrt(k*T/mass);

numP = 100; %number of particles

%box definitions
xmax = 75;
xmin = 0;
ymax = 50;
ymin = 0;
center = xmax/2;


%use 100 steps to get across the region xmax long
t = (200e-9/v_th)/100;


%particle initalization
xmaxNano = 75e-9;
ymaxNano = 50e-9;
 
%  %conductiity equation
%  condEq = charge^2*TauMN*numP/mass;
%  background = 1;
 
 %initalize regions
 [conductivtyN, posN, velN] = regionInit( 'electron', numP);
 [conductivtyP, posP, velP] = regionInit( 'hole', numP);
 

%% Main loops for producing the "movie" of particles


%Probability of scattering
ProbScat = 1- exp(-t/TauMN);
scatterTime = zeros(numP, 1);
ProbMoveWhenHappy = 0.1;
iterations = 20;

figure (4)
plot(posN(:, 1), posN(:, 2), '.b')
hold on
plot(posP(:, 1), posN(:, 2), '.r')
hold off
pause(0.2)
title ('Simulation of Electron Trajectories')
axis([xmin, xmaxNano, ymin, ymaxNano])


for iter =1:iterations
    scatterTime= scatterTime+t*iter;
    
    %Probability of electrons scattering - should it move and how much?
    
    velN = ShouldItMove( t, posN, chargeN, center, ProbScat, ProbMoveWhenHappy, velN, mass, posP, chargeP);
    velP = ShouldItMove( t, posP, chargeP, center, ProbScat, ProbMoveWhenHappy, velP, mass, posN, chargeN);
    
    %Boundary conditions
    
    posN = updatePosition(v_th, numP, posN, velN, t, xmin, xmaxNano, ymin, ymaxNano, chargeN);
    posP = updatePosition(v_th, numP, posP, velP, t, xmin, xmaxNano, ymin, ymaxNano, chargeP);
    
    
    
    figure (4)
    axis([xmin, xmaxNano, ymin, ymaxNano])
    plot(posN(:, 1), posN(:, 2), '.b')
    hold on
    plot(posP(:, 1), posN(:, 2), '.r')
%     hold off
    pause(0.5)
    title ('Simulation of Electron Trajectories')
    
    %Plot updated conductivity fo rmoving electrons in each resion
    
    CondMapUpdated = ConductivityCal(posN, 'electron', numP);
    
    figure(5)
    subplot(1, 2, 1)
    surf(CondMapUpdated)
    title ('Conductivity map changes with time')
    colorbar
    subplot(1, 2, 2)
    plot(posN(1:length(posN), 1), posN(1:length(posN), 2), '.b')
    hold on
    plot(posP(1:length(posP), 1), posN(1:length(posP), 2), '.r')
    axis([xmin, xmaxNano, ymin, ymaxNano])
    pause(0.2)
    title ('Simulation of Electron Trajectories')
    
    
    %temperature map
    %convert velocieites into temperatures, then use hist3 to bin and plot
    Temperature = (velN.^2).*(mass/k);
    
    %     figure(4)
    %     tempDist = hist3(Temperature, [50, 50]);
    %     pcolor(tempDist')
    %     title ('Electron temerature map')
    %     colorbar
    
    
    %calculate the mean free path of the electrons. The time between
    %collisions  is incimetned each iteratin at the top of the iter loop.
    MFP = mean(scatterTime(:, 1));
    
end
