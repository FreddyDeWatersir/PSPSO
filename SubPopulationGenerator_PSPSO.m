function [Swarm, Problem] = SubPopulationGenerator_PSPSO(Dimension, MinCoordinate, MaxCoordinate, InitSwarm, Problem, SwarmSize)
    
    % Evaluate fitness for all particles first
    [fitnessValues, Problem] = fitness(InitSwarm.X, Problem);
    InitSwarm.FitnessValues = fitnessValues;

    % Sort particles by their fitness in descending order
    [~, SortIndex] = sort(fitnessValues, 'descend');

    % Initialize clusters
    clusters = {};
    used = false(size(InitSwarm.X, 1), 1);

    % Form clusters
    for i = 1:size(InitSwarm.X, 1)
        if ~used(SortIndex(i))
            currentCluster = SortIndex(i);
            used(SortIndex(i)) = true;
            distances = pdist2(InitSwarm.X(SortIndex(i), :), InitSwarm.X(~used, :));
            [~, closestIndices] = sort(distances);
            availableIndices = find(~used);
            for j = 1:min(SwarmSize - 1, length(closestIndices))
                currentCluster = [currentCluster; availableIndices(closestIndices(j))];
                used(availableIndices(closestIndices(j))) = true;
            end
            clusters{end+1} = currentCluster;
        end
    end

    % Initialize the swarm structure
    population = struct('X', [], 'Velocity', [], 'FitnessValue', [], 'PbestPosition', [], 'IsConverged', [], 'PbestValue', [], 'GbestValue', [], 'GbestID', [], 'GbestPosition', [], 'Center', [], 'InitRadius', [], 'CurrentRadius', []);
    Swarm = repmat(population, [length(clusters), 1]);

    % Populate the swarm with the clusters
    for i = 1:length(clusters)
        clusterIndices = clusters{i};
        Swarm(i).X = InitSwarm.X(clusterIndices, :);
        Swarm(i).Velocity = -(MaxCoordinate - MinCoordinate) / 4 + (2 * (MaxCoordinate - MinCoordinate) / 4) * rand(size(Swarm(i).X, 1), Dimension);
        %Swarm(i).Shifts = [];
        Swarm(i).FitnessValue = fitnessValues(clusterIndices);
        Swarm(i).PbestPosition = Swarm(i).X;
        Swarm(i).IsConverged = 0;
        Swarm(i).Center = mean(Swarm(i).X, 1);
        Swarm(i).InitRadius = mean(sqrt(sum((Swarm(i).X - Swarm(i).Center).^2, 2)));
        Swarm(i).CurrentRadius = Swarm(i).InitRadius;

        if Problem.RecentChange == 0
            Swarm(i).PbestValue = Swarm(i).FitnessValue;
            [Swarm(i).GbestValue, Swarm(i).GbestID] = max(Swarm(i).PbestValue);
            Swarm(i).GbestPosition = Swarm(i).PbestPosition(Swarm(i).GbestID, :);
        else
            Swarm(i).FitnessValue = -inf(size(Swarm(i).X, 1), 1);
            Swarm(i).PbestValue = Swarm(i).FitnessValue;
            [Swarm(i).GbestValue, Swarm(i).GbestID] = max(Swarm(i).PbestValue);
            Swarm(i).GbestPosition = Swarm(i).PbestPosition(Swarm(i).GbestID, :);
        end
    end
end
