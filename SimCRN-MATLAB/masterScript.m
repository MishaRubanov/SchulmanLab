clear all %#ok<CLALL>
close all
clc

%INPUT VARIABLES
runTime=2*3600;%time to run the simulation for
fontSize=21;
%input rate constants
kf1=0.01;
kr1=0.0025;
kf2=0.001;
kr2=0;
%input initial conditions
a0=1;
b0=1.5;
e0=2;

%maxStep=60;%max step size for ode solver (only needed if ode23s is used

%Construct CRN
crn=simCRN();

%Define Reactions
crn.addRxn({'a','b'},{'c','d'},kf1,kr1);%a+b <-kr1-kf1-> c+d
crn.addRxn({'d','e'},{'f'},kf2,kr2);%d+e <-kr2-kf2-> f
    
%set initial conditions
crn.setConcentration('a',a0);
crn.setConcentration('b',b0);
crn.setConcentration('e',e0);
        
%run the simulation
crn.runSim(runTime,@ode45);
%crn.runSim(runTime,@ode23s,maxStep);%for stiff equations use this

%plot species
species2Plot={'a','b','c','d','e'};
figure('Position', [10 10 6*300 1.25*300]);
hold on
for i=1:length(species2Plot)
    plot(crn.time/3600,crn.conc(crn.getSpeciesIdsByNames(species2Plot{i}),:),'LineWidth',3);
end

z=crn.printCRN();
legend(species2Plot);
xlabel('time (hr)');
ylabel('Concentration (uM)');
set(gca, 'FontSize',fontSize);