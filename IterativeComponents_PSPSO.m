function [Optimizer , Problem] = IterativeComponents_PSPSO(Optimizer,Problem)

    
%% Sub-swarm movements
for ii=1 : Optimizer.SwarmNumber

    if(Optimizer.pop(ii).IsConverged == 0)
        Optimizer.pop(ii).Velocity = Optimizer.w * (Optimizer.pop(ii).Velocity + (Optimizer.c1 * rand(size(Optimizer.pop(ii).X,1) , Optimizer.Dimension).*(Optimizer.pop(ii).PbestPosition - Optimizer.pop(ii).X)) + (Optimizer.c2*rand(size(Optimizer.pop(ii).X,1), Optimizer.Dimension).*(repmat(Optimizer.pop(ii).GbestPosition,size(Optimizer.pop(ii).X,1),1) - Optimizer.pop(ii).X)));
        Optimizer.pop(ii).X = Optimizer.pop(ii).X + Optimizer.pop(ii).Velocity;
        for jj=1 : size(Optimizer.pop(ii).X,1)
            for kk=1 : Optimizer.Dimension
                if Optimizer.pop(ii).X(jj,kk) > Optimizer.MaxCoordinate
                    Optimizer.pop(ii).X(jj,kk) = Optimizer.MaxCoordinate;
                    Optimizer.pop(ii).Velocity(jj,kk) = 0;
                elseif Optimizer.pop(ii).X(jj,kk) < Optimizer.MinCoordinate
                    Optimizer.pop(ii).X(jj,kk) = Optimizer.MinCoordinate;
                    Optimizer.pop(ii).Velocity(jj,kk) = 0;
                end
            end
        end
        [tmp,Problem] = fitness(Optimizer.pop(ii).X,Problem);
        if Problem.RecentChange == 1
            return;
        end
        Optimizer.pop(ii).FitnessValue = tmp;
        for jj=1 : size(Optimizer.pop(ii).X,1)
            if Optimizer.pop(ii).FitnessValue(jj) > Optimizer.pop(ii).PbestValue(jj)
                Optimizer.pop(ii).PbestValue(jj) = Optimizer.pop(ii).FitnessValue(jj);
                Optimizer.pop(ii).PbestPosition(jj,:) = Optimizer.pop(ii).X(jj,:);
                
                [BestPbestValue,BestPbestID] = max(Optimizer.pop(ii).PbestValue);
                if BestPbestValue>Optimizer.pop(ii).GbestValue
                    Optimizer.pop(ii).GbestValue = BestPbestValue;
                    Optimizer.pop(ii).GbestPosition = Optimizer.pop(ii).PbestPosition(BestPbestID,:);
                    Optimizer.pop(ii).GbestID = BestPbestID;
                end
            end
        end
    end
end
%% Update swarm center
for ii=1 : Optimizer.SwarmNumber
    if(Optimizer.pop(ii).IsConverged == 0)
         Optimizer.pop(ii).Center = zeros(1,size(Optimizer.pop(ii).X,2));
         for kk = 1:size(Optimizer.pop(ii).X,2)
            for jj = 1:size(Optimizer.pop(ii).X,1)
                Optimizer.pop(ii).Center(kk) = Optimizer.pop(ii).Center(kk) + Optimizer.pop(ii).PbestPosition(jj,kk);
            end
            Optimizer.pop(ii).Center(kk) = Optimizer.pop(ii).Center(kk)/size(Optimizer.pop(ii).PbestPosition,1);
         end
    end
end
%% Check overlapping and remove worst subpopulation
idx = inf;
while(idx ~= -1)
    idx = -1;
    for ii=1 : Optimizer.SwarmNumber
        if(size(Optimizer.pop(ii).X,1) == 0 || Optimizer.pop(ii).IsConverged == 1) 
           continue;
        end
        for jj=ii+1 : Optimizer.SwarmNumber
            if(size(Optimizer.pop(jj).X,1) == 0 || Optimizer.pop(jj).IsConverged == 1) 
                continue;
            end
            dist = sqrt(sum((Optimizer.pop(ii).GbestPosition - Optimizer.pop(jj).GbestPosition).^2));
            if(dist < Optimizer.pop(ii).InitRadius && dist < Optimizer.pop(jj).InitRadius)
                if(Optimizer.pop(ii).GbestValue > Optimizer.pop(jj).GbestValue)
                    Optimizer.pop(jj) = [];
                    Optimizer.SwarmNumber = Optimizer.SwarmNumber - 1;
                else
                    Optimizer.pop(ii) = [];
                    Optimizer.SwarmNumber = Optimizer.SwarmNumber - 1;
                end
                idx = ii;
                break;
            end
        end
        if(idx ~= -1) 
            break;
        end
    end
end

%% Random Subpop Perturbation
randomSubPopIndex = randperm(Optimizer.SwarmNumber, 1);
[Optimizer.pop(randomSubPopIndex).PbestValue, Problem] = fitness(Optimizer.pop(randomSubPopIndex).PbestPosition, Problem);
[Optimizer.pop(randomSubPopIndex).GbestValue, BestPbestID] = max(Optimizer.pop(randomSubPopIndex).PbestValue);
Optimizer.pop(randomSubPopIndex).GbestPosition = Optimizer.pop(randomSubPopIndex).PbestPosition(BestPbestID, :);
Optimizer.pop(randomSubPopIndex).Velocity = Optimizer.pop(randomSubPopIndex).Velocity + (-Optimizer.PerturbationRange + 2 * Optimizer.PerturbationRange .* rand(size(Optimizer.pop(randomSubPopIndex).X, 1), Optimizer.Dimension));


%% Update Current Radius
for ii=1 : Optimizer.SwarmNumber
    if(Optimizer.pop(ii).IsConverged == 0)
        CurrentRadius = 0.0;
        for jj=1 : size(Optimizer.pop(ii).PbestPosition,1)
            CurrentRadius = CurrentRadius + sqrt(sum((Optimizer.pop(ii).PbestPosition(jj,:) - Optimizer.pop(ii).Center).^2));
        end
        CurrentRadius = CurrentRadius / size(Optimizer.pop(ii).PbestPosition,1);
        Optimizer.pop(ii).CurrentRadius = CurrentRadius;
    end
end

%% Convergence Detection and Deactivation 
AnyConverged = 0;
BestID = GetBestChildID(Optimizer);
for ii=1 : Optimizer.SwarmNumber
    if Optimizer.pop(ii).CurrentRadius < Optimizer.ConvergenceLimit && ii ~= BestID
        Optimizer.pop(ii).IsConverged = 1;
        AnyConverged = AnyConverged + 1;
    end
end

%% Diversity Check and Mechanism
SurvivedParticles = 0;
for ii = 1:Optimizer.SwarmNumber
    if(~Optimizer.pop(ii).IsConverged)
        SurvivedParticles = SurvivedParticles + size(Optimizer.pop(ii).X,1);
    end
end
count = 0;
SaveBestPosition = zeros(AnyConverged,Optimizer.Dimension);
if(SurvivedParticles < Optimizer.initPopulationSize * Optimizer.DiversityDegree)
    while(AnyConverged)
        for ii = 1:Optimizer.SwarmNumber
            if(Optimizer.pop(ii).IsConverged)
                count = count + 1;
                SaveBestPosition(count,:) = Optimizer.pop(ii).GbestPosition;
                Optimizer.pop(ii) = [];
                Optimizer.SwarmNumber = Optimizer.SwarmNumber - 1;
                break;
            end
        end
        AnyConverged = AnyConverged - 1;
    end    

    NumAddParticles = Optimizer.initPopulationSize - SurvivedParticles - count; 
    AddParticles.X = Optimizer.MinCoordinate + ((Optimizer.MaxCoordinate-Optimizer.MinCoordinate).*rand(NumAddParticles,Optimizer.Dimension));
    AddParticles.X = [AddParticles.X;SaveBestPosition];
    [Swarms,Problem] = SubPopulationGenerator_PSPSO(Optimizer.Dimension,Optimizer.MinCoordinate,Optimizer.MaxCoordinate,AddParticles,Problem, Optimizer.SwarmSize);
    if Problem.RecentChange == 1
        return;
    end
    NewSwarms = [Optimizer.pop;Swarms];
    Optimizer.pop = NewSwarms;
    Optimizer.SwarmNumber = length(Optimizer.pop);
end
end

%% Get Best Child ID
function ID = GetBestChildID(Optimizer)
GbestValue = zeros(length(Optimizer.pop),1);
for i = 1:length(Optimizer.pop)
    GbestValue(i) = Optimizer.pop(i).GbestValue;
end
[~,ID] = max(GbestValue);
end