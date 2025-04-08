
function [Problem,E_bbc,E_o,T_r,CurrentError,VisualizationInfo,Iteration] = main_PSPSO(VisualizationOverOptimization,PeakNumber,ChangeFrequency,Dimension,ShiftSeverity,EnvironmentNumber,RunNumber,BenchmarkName)
BestErrorBeforeChange = NaN(1,RunNumber);
OfflineError = NaN(1,RunNumber);
CurrentError = NaN (RunNumber,ChangeFrequency*EnvironmentNumber);
Runtime = NaN(1,RunNumber);
%Algorithm Parameters set up for Excel:

for RunCounter=1 : RunNumber
    if VisualizationOverOptimization ~= 1
        rng(RunCounter+1);%This random seed setting is used to initialize the Problem
    end
    Problem = BenchmarkGenerator(PeakNumber,ChangeFrequency,Dimension,ShiftSeverity,EnvironmentNumber,BenchmarkName);
    rng('shuffle');%Set a random seed for the optimizer
    tic;
    %% Initialiing Optimizer
    clear Optimizer;
    Optimizer.Dimension = Problem.Dimension;
    Optimizer.MaxCoordinate = Problem.MaxCoordinate;
    Optimizer.MinCoordinate = Problem.MinCoordinate;
    Optimizer.w = 0.6;
    Optimizer.c1 = 2.83;
    Optimizer.c2 = 2.83;
    % Optimizer.w = 0.729843788;
    % Optimizer.c1 = 2.05;
    % Optimizer.c2 = 2.05;
    
    %Algorithm Parameters:
    Optimizer.SwarmNumber = 10;
    Optimizer.SwarmSize = 7;
    Optimizer.ConvergenceLimit = 0.01;
    Optimizer.DiversityDegree = 0.7; 
    Optimizer.PerturbationFactor = 0.025; %best 0.025
    Optimizer.ConvergenceLimit = Optimizer.ConvergenceLimit * sqrt(Optimizer.Dimension);
    %Optimizer.OverlapDegree = OptimizerX.OverlapDegree; %originally 0.7
    %Calculate Secondary Parameters
    Optimizer.initPopulationSize = Optimizer.SwarmSize * Optimizer.SwarmNumber;
    Optimizer.PerturbationRange = Optimizer.PerturbationFactor * (Optimizer.MaxCoordinate - Optimizer.MinCoordinate);
    
    InitSwarm.X = Optimizer.MinCoordinate + ((Optimizer.MaxCoordinate-Optimizer.MinCoordinate).*rand(Optimizer.initPopulationSize,Optimizer.Dimension));
    [Optimizer.pop,Problem] = SubPopulationGenerator_PSPSO(Optimizer.Dimension,Optimizer.MinCoordinate,Optimizer.MaxCoordinate,InitSwarm,Problem, Optimizer.SwarmSize);
    Optimizer.SwarmNumber = length(Optimizer.pop);
    VisualizationFlag=0;
    Iteration=0;
    if VisualizationOverOptimization==1
        VisualizationInfo = cell(1,Problem.MaxEvals);
    else
        VisualizationInfo = [];
    end
    %% main loop
    while 1
        Iteration = Iteration + 1;
        %% Visualization for education module
        if (VisualizationOverOptimization==1 && Dimension == 2)
            if VisualizationFlag==0
                VisualizationFlag=1;
                T = Problem.MinCoordinate : ( Problem.MaxCoordinate-Problem.MinCoordinate)/100 :  Problem.MaxCoordinate;
                L=length(T);
                F=zeros(L);
                for i=1:L
                    for j=1:L
                        F(i,j) = EnvironmentVisualization([T(i), T(j)],Problem);
                    end
                end
            end
            VisualizationInfo{Iteration}.T=T;
            VisualizationInfo{Iteration}.F=F;
            VisualizationInfo{Iteration}.Problem.PeakVisibility = Problem.PeakVisibility(Problem.Environmentcounter,:);
            VisualizationInfo{Iteration}.Problem.OptimumID = Problem.OptimumID(Problem.Environmentcounter);
            VisualizationInfo{Iteration}.Problem.PeaksPosition = Problem.PeaksPosition(:,:,Problem.Environmentcounter);
            VisualizationInfo{Iteration}.CurrentEnvironment = Problem.Environmentcounter;
            counter = 0;
            for ii=1 : Optimizer.SwarmNumber
                for jj=1 :size(Optimizer.pop(ii).PbestPosition,1)
                    counter = counter + 1;
                    VisualizationInfo{Iteration}.Individuals(counter,:) = Optimizer.pop(ii).PbestPosition(jj,:);
                end
            end
            VisualizationInfo{Iteration}.IndividualNumber = counter;
            VisualizationInfo{Iteration}.FE = Problem.FE;
        end
        %% Optimization
        [Optimizer,Problem] = IterativeComponents_PSPSO(Optimizer,Problem);
        if Problem.RecentChange == 1%When an environmental change has happened
            Problem.RecentChange = 0;
           
            VisualizationFlag = 0;
            clc; disp(['Run number: ',num2str(RunCounter),'   Environment number: ',num2str(Problem.Environmentcounter)]);
        end
        if  Problem.FE >= Problem.MaxEvals%When termination criteria has been met
            break;
        end
    end
    %% Runtime, Performance Indicator, and Plot Data Gathering
    % This section is common across all EDOAs and must be executed after each run to gather runtime, performance indicators, and plot data.
    elapsedTime = toc;  % Stop the timer for the current run and record runtime
    Runtime(1,RunCounter) = elapsedTime;  % Store the runtime for the current run
    BestErrorBeforeChange(1,RunCounter) = mean(Problem.Ebbc);  % Calculate and store average best error before each environmental change
    OfflineError(1,RunCounter) = mean(Problem.CurrentError);  % Calculate and store the offline error across all function evaluations
    CurrentError(RunCounter,:) = Problem.CurrentError;  % Record current error values for plotting convergence behavior over time
end
%% Output Preparation: Common Across All EDOAs
% This section gathers and summarizes the results of the experiment, including
% performance indicators (E_bbc for Best Error Before Change, E_o for Offline Error),
% runtime statistics (T_r), and any visualization data if enabled.
E_bbc.mean = mean(BestErrorBeforeChange);
E_bbc.median = median(BestErrorBeforeChange);
E_bbc.StdErr = std(BestErrorBeforeChange)/sqrt(RunNumber);
E_bbc.AllResults = BestErrorBeforeChange;
E_o.mean = mean(OfflineError);
E_o.median = median(OfflineError);
E_o.StdErr = std(OfflineError)/sqrt(RunNumber);
E_o.AllResults =OfflineError;
T_r.mean = mean(Runtime);
T_r.median = median(Runtime);
T_r.StdErr = std(Runtime)/sqrt(RunNumber);
T_r.AllResults = Runtime;
if VisualizationOverOptimization==1
    tmp = cell(1, Iteration);
    for ii=1 : Iteration
        tmp{ii} = VisualizationInfo{ii};
    end
    VisualizationInfo = tmp;
end