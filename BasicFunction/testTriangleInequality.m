function violationLevel = testTriangleInequality(distanceMatrix)
    % Input: distanceMatrix is a square matrix of distances between points
    % Output: violationLevel is a value between 0 and 1, where 0 indicates complete violation, and 1 indicates no violation

    % Get the dimensions of the distance matrix
    [n, m] = size(distanceMatrix);

    % Check if the input matrix is square
    if n ~= m
        error('Input matrix must be square');
    end

    % Initialize variables to store the number of satisfied and total inequalities
    satisfiedInequalities = 0;
    totalInequalities = 0;

    % Iterate through all possible combinations of points, taking advantage of symmetry and skipping diagonal elements
    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                % Increment the total number of inequalities tested
                totalInequalities = totalInequalities + 3;

                % Check if the triangle inequality is satisfied for points i, j, and k
                if distanceMatrix(i, j) + distanceMatrix(j, k) >= distanceMatrix(i, k)
                    satisfiedInequalities = satisfiedInequalities + 1;
                end
                if distanceMatrix(i, k) + distanceMatrix(k, j) >= distanceMatrix(i, j)
                    satisfiedInequalities = satisfiedInequalities + 1;
                end
                if distanceMatrix(j, k) + distanceMatrix(k, i) >= distanceMatrix(j, i)
                    satisfiedInequalities = satisfiedInequalities + 1;
                end
            end
        end
    end

    % Calculate the violation level
    violationLevel = 1- satisfiedInequalities / totalInequalities;
end
