# Solving Equation System Using Successive Over Relaxation Method


def SuccessiveOverRelaxation(originMatrix, originVectorB, W, Accuracy):
    """
    Solving Equation System in the Successive Over Relaxation method

    :param originMatrix: NxN Matrix
    :param originVectorB: Nx1 Vector
    :param W: The relaxation factor
    :param Accuracy: The Accuracy for the solution
    """
    # In case the input equation system isn't Quadratic or in the appropriate size
    if len(originMatrix) != len(originMatrix[0]) or len(originVectorB) != len(originMatrix) or len(originVectorB[0]) != 1:
        printIntoFile(None, "The input equation system isn't match", False)
        print("The input equation system isn't match")

    # In case the matrix has more or less than one solution
    if determinantMatrix(originMatrix) == 0:
        printIntoFile(None, 'This is Singular matrix', False)
        print('This is Singular matrix')

    # In case the omega parameter is out of boundaries
    if W <= 0 or W >= 2:
        printIntoFile(None, 'Omega parameter is out of boundaries', False)
        print('Omega parameter is out of boundaries')

    # Organize the matrix pivots
    originMatrix, originVectorB = organizeMatrix(originMatrix, originVectorB)

    # Our lists for the Preview iteration values, and our Current iteration values
    prevIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]
    currentIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]

    # Loop for finding the solution
    for _ in range(500):

        # Calculate the next guess
        for i in range(len(originMatrix)):
            rowSum = 0
            for j in range(len(originMatrix)):
                if i != j:
                    rowSum = rowSum + originMatrix[i][j] * currentIteration[j][0]

            # Update the new value of x for this iteration
            currentIteration[i][0] = (1 - W) * prevIteration[i][0] + (W * (originVectorB[i][0] - rowSum) / originMatrix[i][i])

        # Save the current iteration values into the file
        printIntoFile(currentIteration, _ + 1, True)

        # In case we found the solution, Stop the program
        if all([False if abs(currentIteration[row][0] - prevIteration[row][0]) > Accuracy else True for row in range(len(currentIteration))]):
            break

        # Update the previous solution to be the current solution
        prevIteration = [[currentIteration[row][0] for _ in range(1)] for row in range(len(currentIteration))]

        # According message in case the Matrix is not converge
        if _ == 499:
            if not isDiagonalDominant(originMatrix):
                printIntoFile(None, "This isn't Diagonal Dominant matrix", False)
                print("This isn't Diagonal Dominant matrix")

            printIntoFile(None, "This equation system isn't Converge", False)
            print("This equation system isn't Converge")
            exit()

    # Saving the equation system final solution
    printIntoFile(currentIteration, 'Solution', True)
    print(f'Equation system solution {list(map(lambda x: int(x[0] * 10 ** solutionPrecision(Accuracy)) / 10 ** solutionPrecision(Accuracy), currentIteration))}')


def organizeMatrix(originMatrix, originVectorB):
    """
    Taking care that the pivot in the every row will be the highest possible, and return the updated equation system

    :param originMatrix: NxN matrix
    :param originVectorB: Nx1 vector
    :return: The updated equation system
    """
    # Saving the equation system the user gave
    EquationSystem = [[originMatrix[row][col] for col in range(len(originMatrix[0]))] for row in range(len(originMatrix))]
    [EquationSystem[row].append(originVectorB[row][0]) for row in range(len(originVectorB))]
    printIntoFile(EquationSystem, 'Inserted Equation System\n', False)

    # Loop to get the highest pivots possible
    for i in range(len(originMatrix)):

        # Variable to store the highest value for the pivot
        maxPivot = abs(originMatrix[i][i])

        # Variable to store the new pivot row
        pivotRow = -1

        # Searching the highest potential Pivot for originMatrix[i][i]
        for j in range(i + 1, len(originMatrix)):

            # In case there's a higher pivot (on the Column[i])
            if abs(originMatrix[j][i]) > maxPivot:
                maxPivot = abs(originMatrix[j][i])
                pivotRow = j

        # In case there was a higher pivot, change the matrix so the Pivot will be the maximum
        if maxPivot != abs(originMatrix[i][i]):
            originVectorB[i], originVectorB[pivotRow] = originVectorB[pivotRow], originVectorB[i]
            originMatrix[i], originMatrix[pivotRow] = originMatrix[pivotRow], originMatrix[i]

    # Saving the equation system after changing rows/cols
    EquationSystem = [[originMatrix[row][col] for col in range(len(originMatrix[0]))] for row in range(len(originMatrix))]
    [EquationSystem[row].append(originVectorB[row][0]) for row in range(len(originVectorB))]
    printIntoFile(EquationSystem, 'Updated Equation System', False)

    # Return the updated equation system
    return originMatrix, originVectorB


def isDiagonalDominant(matrix):
    """
    Check if the pivot in every row is bigger than the sum of the whole row (without the pivot),
    If yes return True, else False

    """
    for i in range(len(matrix)):

        # Variable to store, the summation of absolute row [i]
        rowSum = 0
        for j in range(len(matrix)):
            if i != j:
                rowSum = rowSum + abs(matrix[i][j])

        # If the summation of the row is bigger than the pivot, return False (The matrix is not diagonal dominant)
        if rowSum > abs(matrix[i][i]):
            return False

    # The matrix is Diagonal Dominant
    return True


def determinantMatrix(matrix):
    """
    Calculate the matrix determinant and return the result

    :param matrix: NxN Matrix
    :return: Matrix determinant
    """
    # Simple case, The matrix size is 2x2
    if len(matrix) == 2:
        value = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]
        return value

    # Initialize our sum variable
    determinantSum = 0

    # Loop to traverse each column of the matrix
    for current_column in range(len(matrix)):
        sign = (-1) ** current_column

        # Calling the function recursively to get determinant value of sub matrix obtained
        determinant_sub = determinantMatrix(
            [row[: current_column] + row[current_column + 1:] for row in (matrix[: 0] + matrix[0 + 1:])])

        # Adding the calculated determinant value of particular column matrix to total the determinantSum
        determinantSum = determinantSum + (sign * matrix[0][current_column] * determinant_sub)

    # Returning the final Sum
    return determinantSum


def printIntoFile(data, message, isVector):
    """
    Printing the content into a specified file

    :param data: Data is a list representing matrix
    :param message: Message is a string representing a message
    :param isVector: isVector is a boolean representing if the data is a vector
    """
    # Open file and save the sent content
    with open('Calculation.txt', 'a+') as file:

        # In case we sent a message
        if message:
            file.write('\n{: ^25}'.format(message))
            file.write('' if message != 'Updated Equation System' else '\n')

        # In case we sent a data
        if data:
            for i in range(len(data)):
                for j in range(len(data[i])):
                    file.write('{: ^25}'.format(int(data[i][j] * 10 ** solutionPrecision(Precision)) / 10 ** solutionPrecision(Precision)))
                file.write('' if isVector else '\n')

        # Used to enhance the appearance
        if message == 'Updated Equation System':
            file.write('\n==========================================================================================\n')
            for i in range(len(data) + 1):
                file.write('{: ^25}'.format('Iteration' if i == 0 else chr(64 + i)))


def resetFile():
    """
    Reset the calculation fileS

    """
    with open('Calculation.txt', 'w') as file:
        file.write('------------------------------ Successive Over Relaxation Method ------------------------------')


def solutionPrecision(accuracy):
    """
    Method for getting the right precision of the solution

    :param accuracy: Precision of the solution
    """
    counterPower = 0
    while accuracy != 1:
        accuracy = accuracy * 10
        counterPower = counterPower + 1

    # return the precision after the dot of the solution
    return counterPower


# Our Program Driver
if __name__ == "__main__":

    # Reset the calculation file
    resetFile()

    # Input section
    inputMatrix = [[4, 3, 0], [3, 4, -1], [0, -1, 4]]
    inputVectorB = [[24], [30], [-24]]
    Omega = 1.2
    Precision = 0.0001

    # Running the program
    print('---------- Successive Over Relaxation Method ----------')
    SuccessiveOverRelaxation(inputMatrix, inputVectorB, Omega, Precision)
    print('\n\nCalculation Is Done, Check File "Calculation" For More Information')

