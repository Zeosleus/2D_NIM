# 2-D NIM Game
import os
import random


############################### FG COLOR DEFINITIONS ###############################
class bcolors:
    # pure colors...
    GREY = '\033[90m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BLUE = '\033[94m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    # color styles...
    HEADER = '\033[95m'
    QUESTION = '\033[93m\033[3m'
    MSG = '\033[96m'
    WARNING = '\033[93m'
    ERROR = '\033[91m'
    ENDC = '\033[0m'  # RECOVERS DEFAULT TEXT COLOR
    BOLD = '\033[1m'
    ITALICS = '\033[3m'
    UNDERLINE = '\033[4m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


def screen_clear():
    # for mac and linux(here, os.name is 'posix')
    if os.name == 'posix':
        _ = os.system('clear')
    else:
        # for windows platform
        _ = os.system('cls')


def initializeBoard(N):
    board = [''] * (N * N + 1)

    # this is the COUNTER of cells in the board already filled with R or G
    board[0] = 0

    # each EMPTY cell in the board contains its cardinal number
    for i in range(N * N):
        if i < 9:
            board[i + 1] = ' ' + str(i + 1)
        else:
            board[i + 1] = str(i + 1)

    return board


def drawNimPalette(board, N):
    EQLINE = '\t'
    MINUSLINE = '\t'
    CONSECEQUALS = ''
    CONSECMINUS = ''
    for i in range(5):
        CONSECEQUALS = CONSECEQUALS + '='
        CONSECMINUS = CONSECMINUS + '-'

    for i in range(10):
        EQLINE = EQLINE + CONSECEQUALS
        MINUSLINE = MINUSLINE + CONSECMINUS

    for i in range(N):
        # PRINTING ROW i...
        if i == 0:
            print(EQLINE)
        else:
            print(MINUSLINE)

        printRowString = ''

        for j in range(N):
            # PRINTING CELL (i,j)...
            CellString = str(board[N * i + j + 1])
            if CellString == 'R':
                CellString = ' ' + bcolors.RED + CellString + bcolors.ENDC

            if CellString == 'G':
                CellString = ' ' + bcolors.GREEN + CellString + bcolors.ENDC

            if printRowString == '':
                printRowString = '\t[ ' + CellString
            else:
                printRowString = printRowString + ' | ' + CellString
        printRowString = printRowString + ' ]'
        print(printRowString)
    print(EQLINE)
    print(bcolors.PURPLE + '\t\t\tCOUNTER = [ ' + str(board[0]) + ' ]' + bcolors.ENDC)
    print(EQLINE)


def inputPlayerLetter():
    # The player chooses which label (letter) will fill the cells
    letter = ''
    while not (letter == 'G' or letter == 'R'):
        print(bcolors.QUESTION + '[Q1] What letter do you choose to play? [ G(reen) | R(ed) ]' + bcolors.ENDC)
        letter = input().upper()
        # The first letter corresponds to the HUMAN and the second element corresponds to the COMPUTER
        if letter == 'G':
            return ['G', 'R']
        else:
            if letter == 'R':
                return ['R', 'G']
            else:
                print(bcolors.ERROR + 'ERROR1: You provided an invalid choice. Please try again...' + bcolors.ENDC)


def whoGoesFirst():
    if random.randint(0, 1) == 0:
        return 'computer'
    else:
        return 'player'


def howComputerPlays():
    while True:
        print(
            bcolors.QUESTION + '[Q5] How will the computer play? [ R (randomly) | F (first Free) | C (copycat)]' + bcolors.ENDC)
        strategyLetter = input().upper()

        if strategyLetter == 'R':
            return 'random'
        else:
            if strategyLetter == 'F':
                return 'first free'
            else:
                if strategyLetter == 'C':
                    return 'copycat'
                else:
                    print(
                        bcolors.ERROR + 'ERROR 3: Incomprehensible strategy was provided. Try again...' + bcolors.ENDC)


def getBoardSize():
    BoardSize = 0
    while BoardSize < 1 or BoardSize > 10:
        GameSizeString = input('Determine the size 1 =< N =< 10, for the NxN board to play: ')
        if GameSizeString.isdigit():
            BoardSize = int(GameSizeString)
            if BoardSize < 1 or BoardSize > 10:
                print(
                    bcolors.ERROR + 'ERROR 4: Only positive integers between 1 and 10 are allowable values for N. Try again...' + bcolors.ENDC)
        else:
            print(
                bcolors.ERROR + 'ERROR 5: Only positive integers between 1 and 10 are allowable values for N. Try again...' + bcolors.ENDC)
    return (BoardSize)


def startNewGame():
    # Function for starting a new game
    print(bcolors.QUESTION + '[Q0] Would you like to start a new game? (yes or no)' + bcolors.ENDC)
    return input().lower().startswith('y')


def continuePlayingGame():
    # Function for starting a new game
    print(bcolors.QUESTION + '[Q2] Would you like to continue playing this game? (yes or no)' + bcolors.ENDC)
    return input().lower().startswith('y')


def playAgain():
    # Function for replay (when the player wants to play again)
    print(bcolors.QUESTION + '[Q3] Would you like to continue playing this game? (yes or no)' + bcolors.ENDC)
    return input().lower().startswith('y')


def isBoardFull(board, N):
    return board[0] == N * N


# Helper function to check if the current board is symmetric.
# A board is symmetric iff for every empty cell (i, j), its symmetric cell (j, i) is also empty.
def isBoardSymmetric(nimBoard, emptyCells, N):
    emptySymmetricCellCount = 0

    for cell in emptyCells:
        cellRow, cellCol = getRowAndColumn(cell, N)
        if isCellDiagonal(cell, N):
            symmetricCellIndex = (N - cellRow) * N + (N - cellCol + 1)
        else:
            symmetricCellIndex = (cellCol - 1) * N + cellRow

        if isCellEmpty(nimBoard, cell) and isCellEmpty(nimBoard, symmetricCellIndex):
            emptySymmetricCellCount = emptySymmetricCellCount + 1

    return True if emptySymmetricCellCount == len(emptyCells) else False


# Helper function to check if the given cell is empty or not.
def isCellEmpty(nimBoard, cellIndex):
    return not (nimBoard[cellIndex] == 'G' or nimBoard[cellIndex] == 'R')


# Helper function to check if a given cell is a diagonal cell.
# The indexes of the diagonal cells, differ by N + 1, so cell k is diagonal, iff : (k - 1) mod (N + 1) = 0
def isCellDiagonal(cellIndex, N):
    return (cellIndex - 1) % (N + 1) == 0


def getRowAndColumn(move, N):  # returns the row and col of the "move" argument
    moveRow = 1 + (move - 1) // N
    moveColumn = 1 + (move - 1) % N
    return moveRow, moveColumn


# Helper function : returns a list with the linear indexes of the diagonal cells of the NxN board
def getDiagCells(N):
    diagCells = []
    index = 1
    while index <= N * N:
        diagCells.append(index)
        index = index + N + 1

    return diagCells


# Helper function to return random number of additional cells belonging to the same row or col as firstCell
# and such that all cells when put together, form a block of consecutive cells.
def findRandomCandidateCells(nimBoard, firstCell):
    rndDir = random.choice(['row', 'col'])  # select either row or column direction
    rndNumCells = random.randint(0, maxNumMoves - 1)  # select the random number of additional cells

    selectedCells = [firstCell]  # the final selected cells that make-up the random move
    candidateCells = [firstCell]  # the candidate cells with respect to the first random cell (firstCell) and the rndDir
    firstCellRow, firstCellCol = getRowAndColumn(firstCell, N)

    if rndDir == "row":
        loopStartIndex = (firstCellRow - 1) * N + 1  # the linear index of the first cell in firstCell's row
        loopEndIndex = loopStartIndex + N  # the linear index of the last cell in firstCell's row
        loopStep = 1  # run-through a board line
        distThresh = rndNumCells  # we will keep the rndNumCells left and right cells
        consecutiveCellDist = 1  # in "row" mode, the linear indexes of consecutive cells, differ by 1 (abs value)
    elif rndDir == "col":
        loopStartIndex = firstCellCol  # the linear index of the first cell in firstCell's col
        loopEndIndex = loopStartIndex + N * N  # the linear index of the last cell in firstCell's col
        loopStep = N  # run-through a board column
        distThresh = rndNumCells * N  # we will keep the rndNumCells previous and next cells, of the same column
        consecutiveCellDist = N  # in "col" mode, the linear indexes of consecutive cells, differ by N (abs value)

    # Gather all the other (i != firstCell, firstCell already belongs to the candidateCells list), non-diagonal and
    # empty cells, that belong to the same row/col as firstCell.
    # The candidate cells must be at a distance of "rndNumCells cells" to the left/right or above/below.
    for i in range(loopStartIndex, loopEndIndex, loopStep):
        if i != firstCell and not isCellDiagonal(i, N) and isCellEmpty(nimBoard, i) \
                and abs(firstCell - i) <= distThresh:
            candidateCells.append(i)

    # Sort the candidateCells list, in order to only keep the cells that form a consecutive block, when examined from
    # left (cell with the smallest linear index) to right (cell with the largest linear index)
    candidateCells.sort()

    # If at least one candidate cell was found
    if len(candidateCells) > 1:
        # Remove cells whose distance from the previous cell, is greater than the consecutive cell distance
        # (i.e. 1 for "row" mode and N for "col" mode).
        # The reason that firstCell must be a part of the candidateCells list, is to perform the following check,
        # regarding the distance of the consecutive cells.
        for i in range(1, len(candidateCells)):
            if abs(candidateCells[i] - candidateCells[i - 1]) > consecutiveCellDist:
                candidateCells.pop(i)

        # Remember that candidateCells contains firstCell, so after the above "non-consecutive cell elimination" process
        # we need to remove it
        candidateCells.remove(firstCell)
        random.shuffle(candidateCells)  # permute the candidateCells list
        for i in range(0, rndNumCells):  # pick a random candidate cell
            if len(candidateCells) > 0:
                # remove the last element from the candidateCells list (in order to avoid duplicate random cells)
                # and append it to the selectedCells list
                selectedCells.append(candidateCells.pop())

    return selectedCells


# Implements the 'Random' computer strategy
def getComputerMove_random(nimBoard, emptyCells, N):
    # The computer must randomly choose a first empty cell.
    # If that cell is a diagonal one, stop (i.e. choose no other cells).
    # If the chosen cell is NOT a diagonal one, randomly choose between rows and columns.
    # Then, choose 0 to 2 more cells, such that:
    #     * They belong to the same row/col with the first randomly chosen cell
    #     * They are not diagonal cells
    #     * All the cells are consecutive
    firstCell = random.choice(emptyCells)  # choose a random **empty** cell
    chosenCells = [firstCell]

    # If the random cell is non-diagonal, select a random number of additional cells (0 to 2), either in row or
    # column direction
    if not isCellDiagonal(firstCell, N):
        # The list returned by the following function, contains firstCell as well as the other randomly selected
        # additional cells
        chosenCells = findRandomCandidateCells(nimBoard, firstCell)

    return chosenCells


# Implements the 'Firstfit' computer strategy
def getComputerMove_firstfit(nimBoard, emptyCells, N):
    # The computer chooses as the move's first cell, the first empty cell.
    # If that cell is a diagonal one, stop (i.e. choose no other cells).
    # If the chosen cell is NOT a diagonal one, randomly choose between rows and columns.
    # Then, choose 0 to 2 more cells, such that:
    #     * They belong to the same row/col with the first randomly chosen cell
    #     * They are not diagonal cells
    #     * All the cells are consecutive
    emptyCells.sort()  # sort the emptyCells list -> emptyCells[0] = index of the first empty board cell
    firstEmptyCell = emptyCells[0]  # get the index of the first empty cell
    chosenCells = [firstEmptyCell]

    # If the first empty cell found is non-diagonal, select a random number of additional cells (0 to 2), either in row
    # or column direction
    if not isCellDiagonal(firstEmptyCell, N):
        # The list returned by the following function, contains firstCell as well as the other randomly selected
        # additional cells
        chosenCells = findRandomCandidateCells(nimBoard, firstEmptyCell)

    return chosenCells


# Implements the 'Copycat' computer strategy
def getComputerMove_copycat(nimBoard, emptyCells, opponentMove, N):
    chosenCells = []

    # If the computer get chosen to play at the very first round, there is no prior HUMAN move to copy, so we
    # assume that the first move by the computer is made using the 'firstfit' strategy.
    if len(opponentMove) == 0:
        chosenCells = getComputerMove_firstfit(nimBoard, emptyCells, N)

    else:
        opponentMoveCellCoords = []
        # Each element of the list will be a tuple with the 2D coordinates of each cell in the HUMAN player's last move
        for cell in opponentMove:
            opponentMoveCellCoords.append(getRowAndColumn(cell, N))

        # If the first cell in the HUMAN player's last move, is a *diagonal* one, that means that the move only
        # consisted of a diagonal cell, since if a diagonal cell is chosen,
        # no others cells can be player during the same move.
        if isCellDiagonal(opponentMove[0], N):  # try to play the next diagonal cell
            nextDiagCellRow = N - opponentMoveCellCoords[0][0] + 1
            nextDiagCellCol = N - opponentMoveCellCoords[0][1] + 1
            nextDiagCell = (nextDiagCellRow - 1) * N + nextDiagCellCol
            if isCellEmpty(nimBoard, nextDiagCell):  # play that diagonal cell
                chosenCells.append(nextDiagCell)
            else:
                diagCells = getDiagCells(N)  # a list with the linear indexes of all the diagonal cells
                occupiedCount = 0
                for cell in diagCells:
                    if not isCellEmpty(nimBoard, cell):
                        occupiedCount = occupiedCount + 1

                if occupiedCount == len(diagCells):  # check if all the diagonal cells are occupied
                    emptyCells.sort()
                    chosenCells.append(emptyCells[0])  # all diagonal cells occupied -> chose the first empty cell
                elif occupiedCount < len(diagCells):
                    diagCells.remove(nextDiagCell) # remove the desired but occupied diagonal cell
                    # choose a random *empty diagonal* cell to play
                    random.shuffle(diagCells)  # randomly permute the diagCells list
                    while True:
                        rndIndex = random.randint(0, len(diagCells) - 1)  # random index
                        if isCellEmpty(nimBoard, diagCells[rndIndex]):
                            chosenCells.append(diagCells[rndIndex])
                            break
                        else:
                            diagCells.remove(rndIndex)

        else:
            # Each element of the opponentMoveCellCoords list is a tuple that contains the 2D board coordinates
            # (i.e. (i, j) of each cell in the HUMAN player's last move (opponentMoveCellCoords[k][0] stores the
            # row coordinate ('i') and opponentMoveCellCoords[k][1] stores the column coordinate ('j') for cell
            # number k, in the HUMAN player's last move).
            # In the copycat strategy the COMPUTER must play the cells that are symmetric to the cells the HUMAN player
            # chose, with respect to the board's main diagonal.
            # Therefore, if the HUMAN player chose the cell (i, j), the COMPUTER must choose the cell (j, i).
            # Remember that cell (i, j) corresponds to the linear index (i - 1) * N + j.
            # Thus, for each move the HUMAN player made, we reverse the coordinates from (i, j) to (j, i), and then we
            # compute the linear index that corresponds to the 'copycat' cell (j, i).
            for i in range(0, len(opponentMoveCellCoords)):
                chosenCells.append((opponentMoveCellCoords[i][1] - 1) * N + opponentMoveCellCoords[i][0])

    return chosenCells


# Helper function to check if the cells in the 'cells' argument, form a block of consecutive cells
# either in row or col direction.
def areCellsConsecutive(cells, N):
    sameRowCount = 0
    sameColCount = 0
    successiveRowCount = 0
    successiveColCount = 0

    # If we don't sort the input cells list, if, for example on a 4x4 board, cells 13, 14, 15 are all empty,
    # when the user enters 13 as the first cell and then 15 14, this will be treated as an illegal move,
    # even if it is actually allowed, because the following code segment that checks if all cells form a
    # consecutive block, compares the "distance" of each cell from its previous one.
    # Therefore, in the example mentioned above, the input list would be [13, 15, 14] and abs(15-13) = 2 > 1.
    # So, if we sort the list prior to the "all cells consecutive" check, we account for the case where
    # the user enters 13, 15, 14 as well as the "sorted" case of 13, 14, 15.
    cells.sort()

    # Check if the i-th cell (i >= 1) belongs to the same row/col as cells[0], i.e. the first cell.
    # The relation "in same row/col" is transitive, therefore if cells[0] and cells[1] belong to the same row/col
    # and cells[0] and cells[2] belong to the same row/col, then cells[1] and cells[2] also belong to the same row/col.
    # Also, check if neither cell is diagonal.
    for i in range(1, len(cells)):
        if getRowAndColumn(cells[0], N)[0] == getRowAndColumn(cells[i], N)[0] \
                and not isCellDiagonal(cells[0], N) and not isCellDiagonal(cells[i], N):
            sameRowCount = sameRowCount + 1
        elif getRowAndColumn(cells[0], N)[1] == getRowAndColumn(cells[i], N)[1] \
                and not isCellDiagonal(cells[0], N) and not isCellDiagonal(cells[i], N):
            sameColCount = sameColCount + 1

    if sameRowCount == len(cells) - 1 or sameColCount == len(cells) - 1:
        # check if cells form a consecutive block
        for i in range(0, len(cells) - 1):
            if abs(cells[i] - cells[i + 1]) == 1:  # the indexes of successive cells of the same row, differ by 1
                successiveRowCount = successiveRowCount + 1
            elif abs(cells[i] - cells[i + 1]) == N:  # the indexes of successive cells of the same col, differ by N
                successiveColCount = successiveColCount + 1

    if (sameRowCount == len(cells) - 1 or sameColCount == len(cells) - 1) \
            and (successiveRowCount == len(cells) - 1 or successiveColCount == len(cells) - 1):
        return True
    else:
        return False


# Helper function to return all possible consecutive cell blocks of size blockDim
def getAllPossibleBlocks(emptyCells, blockDim, N):
    blockMoves = []

    # Sort the emptyCells list, since the block detection logic depends on starting from "block-start" cells and
    # examining the next blockDim - 1 "places" (i.e. cells whose linear index is bigger than the "block-start" cell)
    emptyCells.sort()

    # choose "block-start" cell
    for i in range(0, len(emptyCells) - (blockDim - 1)):
        # skip diagonal cells, since no cell block can start from a diagonal cell
        if not isCellDiagonal(emptyCells[i], N):
            # gather all cells in the next blockDim - 1 "places" of the same **row** as the one the "block-start" cell
            # belongs to. All cells must be non-diagonal, since no cell block can contain diagonal cells.
            rowNextCells = [cell for cell in [*range(emptyCells[i], emptyCells[i] + blockDim, 1)]
                            if not isCellDiagonal(cell, N)
                            and getRowAndColumn(emptyCells[i], N)[0] == getRowAndColumn(cell, N)[0]]

            # Since rowNextCells contains emptyCells[i], i.e. the current block-start cell, we need to check if
            # the number of cells in rowNextCells is equal to the number of cells in the block (i.e. blockDim), and
            # also that each cell in rowNextCells is empty.
            # The same applies to colNextCells below.
            if len(rowNextCells) == blockDim and all(cell in emptyCells for cell in rowNextCells):
                blockMoves.append(tuple(emptyCells[i: i + blockDim]))

            # gather all cells in the next blockDim - 1 "places" of the same **column** as the one the "block-start"
            # cell belongs to. All cells must be non-diagonal, since no cell block can contain diagonal cells.
            colNextCells = [cell for cell in [*range(emptyCells[i], emptyCells[i] + blockDim * N, N)]
                            if not isCellDiagonal(cell, N)
                            and getRowAndColumn(emptyCells[i], N)[1] == getRowAndColumn(cell, N)[1]]

            if len(colNextCells) == blockDim and all(cell in emptyCells for cell in colNextCells):
                blockMoves.append(tuple(colNextCells))

    return blockMoves  # return the list that contains both all possible row and col blocks, as a list of "block" moves


# This function is called from getComputerWinningMove, to check if the current board is a *winning* state,
# ONLY IF THERE ARE AT LEAST 3 EMPTY CELLS.
# A winning state is a state where there exists **at least one** of its alternative games (options),
# where that game is a winning game.
# Otherwise, if all possible moves, lead to (simpler) losing states, then the current state is a *losing* state.
# So given the current state of the board -> check all possible alternatives games if they are winning games.
def checkIfWinState(emptyCells, N):
    # Base winning state case of 2 empty and consecutive cells
    if len(emptyCells) == 2 and areCellsConsecutive(emptyCells, N):
        return True, []
    # Base losing state case of 2 empty and non-consecutive cells
    elif len(emptyCells) == 2 and not areCellsConsecutive(emptyCells, N):
        return False, []

    # if the board is symmetric, we are in a losing state (G + G = *0 (copycat principle)),
    # i.e. no need to check the available player options
    if isBoardSymmetric(nimBoard, emptyCells, N):
        return False, []
    else:
        # Try to play 1 cell among the currently empty ones
        for i in range(0, len(emptyCells)):
            currMove = [emptyCells[i]]
            # Assume that the i-th empty cell (emptyCells[i]), was chosen as the current move.
            # Try to see if that move gives a winning state, by checking the alternative games it leads to
            newEmptyCells = [cell for cell in emptyCells if cell != emptyCells[i]]

            # If the previous recursive call returned True that means the current move resulted in a game where the
            # only empty cells were two consecutive ones, which is a base winning state case.
            # Therefore, if after applying the current move, the recursive call returned False, that means that the
            # current move leads to a losing state for the HUMAN player, i.e. a winning state for the COMPUTER.
            if not checkIfWinState(newEmptyCells, N)[0]:
                return True, currMove

        # Try to play 2 cells among the currently empty ones
        # The 2 cells to check are not chosen randomly, they need to form a consecutive block in row or col dir
        # Thus, we must try to find 2 consecutive cells in the row dir or col dir -> call getAllPossibleBlocks()
        # to get all possible blocks of size 2, either in row or col direction, and try to play each one of them.
        cellBlocksSize2 = getAllPossibleBlocks(emptyCells, maxNumMoves - 1, N)
        for block in cellBlocksSize2:
            currMove = block  # try to play this cell block
            newEmptyCells = [cell for cell in emptyCells if cell not in currMove]

            # If the current move leads to a losing state for the opponent -> COMPUTER wins
            if not checkIfWinState(newEmptyCells, N)[0]:
                return True, currMove

        # try to play 3 cells among the currently empty ones
        # The 3 cells to check are not chosen at random, they need to form a consecutive block either in row or col dir.
        # Thus, we must try to find 3 consecutive cells in the row or col dir -> call getAllPossibleBlocks()
        # to get all possible blocks of size 3, either in row or col direction, and try to play each one of them.
        cellBlocksSize3 = getAllPossibleBlocks(emptyCells, maxNumMoves, N)
        for block in cellBlocksSize3:  # try to play the current block
            currMove = block  # try to play this cell block
            newEmptyCells = [cell for cell in emptyCells if cell not in currMove]

            # If the current move leads to a losing state for the opponent -> COMPUTER wins
            if not checkIfWinState(newEmptyCells, N)[0]:
                return True, currMove

        return False, []  # no winning moves were found, i.e. all possible next games correspond to losing states


# Function that returns the optimal move of the COMPUTER player.
# The case where there is only **1** empty cell left, corresponds to a trivial winning state, therefore there is no
# need to call the recursive 'checkIfWinState' function.
# We assume that the case where there are **2** empty cells left, is the base of the induction implemented using the
# recursive 'checkIfWinState' function.
# Specifically, if there are 2 empty and consecutive cells (where neither is a diagonal cell), this is a winning state.
# Otherwise, if there are 2 empty and non-consecutive cells, this corresponds to a losing state.
# For the above cases, there is no need to call the recursive function.
# However, when there are at least **3** empty cells, we use the recursive function to calculate the optimal move for
# the COMPUTER player.
def getComputerWinningMove(emptyCells, N):
    if len(emptyCells) == 1:  # if there is only one empty cell, no point in calling the recursive function
        return True, [emptyCells[0]]
    elif len(emptyCells) == 2:
        if areCellsConsecutive(emptyCells, N):  # two successive cells -> take them and win
            return True, [emptyCells[0]] + [emptyCells[1]]
        else:
            return False, []
    elif len(emptyCells) > 2:  # if there are at least 3 cells, call the recursive function to find the winning move
        return checkIfWinState(emptyCells, N)


def writeCells(nimBoard, emptyCells, cells, letter):
    for cell in cells:
        if isCellEmpty(nimBoard, cell):
            nimBoard[cell] = letter  # "color" the cell using the corresponding player's letter (R/G)
            emptyCells.remove(cell)  # remove cell from the emptyCells list
            nimBoard[0] = nimBoard[0] + 1  # update board counter


# Helper function to check if the cells given by the user, comply with the game's rules
def validateUserInput(nimBoard, inputCells, N):
    # necessary checks :
    # 1) are the given cells empty ? (done inside getAndCheckInputCells())
    # 2) are the given cells non-diagonal
    # 3) are all the cells in the same line/col as the firstCell?
    # 4) are all the cells consecutive
    emptyCount = 0
    nonDiagCellCount = 0
    numNewCells = len(inputCells) - 1

    if numNewCells > 0:  # if the user entered at least one additional cell
        # check 1 : are all cells empty?
        for i in range(1, numNewCells + 1):
            if isCellEmpty(nimBoard, inputCells[i]):
                emptyCount = emptyCount + 1

        # check 2 : are all the cells non-diagonal?
        for i in range(1, numNewCells + 1):
            if not isCellDiagonal(inputCells[i], N):
                nonDiagCellCount = nonDiagCellCount + 1

        # check 3 : are all the cells in the same line/col as the firstCell?
        #                       and
        # check 4 : are all the cells consecutive, are performed inside "areCellsConsecutive()".

        # all checks passed?
        if emptyCount == numNewCells and nonDiagCellCount == numNewCells and areCellsConsecutive(inputCells, N):
            return True
        else:
            print(bcolors.GREY + "The cells you entered, violate the game's rules... Please try again, or press enter "
                                 "to continue using only the first cell entered" + bcolors.ENDC)
            return False

    elif numNewCells == 0:
        return True


def getAndCheckInputCells(nimBoard, mode, cells, maxNumMoves, N):
    if mode == "first":  # the move's first cell
        while True:  # while the first cell is not an empty one, ask for a new input cell
            userInput = input("Enter the linear index of the first **empty** cell : ")
            if userInput != "":  # "enter" not pressed
                if len(userInput) == 1:  # if the user enters more than 2 numbers, "ignore" them
                    firstCell = int(userInput)
                    if isCellEmpty(nimBoard, firstCell):  # check if the given cell is empty
                        break

        return firstCell

    elif mode == "next":  # the move's additional cells -> check their validity using "validateUserInput()".
        while True:
            # Split the input string using whitespace as a delimiter, and then used List Comprehension to construct
            # a list where each element is the integer representation of the number given as a string by the user.
            nextCells = [int(x) for x in
                         input("You may give up to " + str(maxNumMoves - 1) + " two more cells : ").split()]

            # if at most maxNumMoves-1 (i.e. 2 in the default case) additional cells were entered
            if len(nextCells) <= maxNumMoves - 1:
                # When the function is called using the "next" mode, the 'cells' argument contains the first cell
                # given by the user. The 'newCells' list contains the first cell given by the user, as well all the
                # additional cells. This is needed in order for 'validateUserInput()' to be able to check if the
                # additional cells are in the same row or col with the first cell and that all cells form a
                # consecutive block.
                # When 'validateUserInput()' returns, we need to remove the first cell from the list to only return the
                # move's additional cells.
                newCells = cells + nextCells

                if validateUserInput(nimBoard, newCells, N):
                    break
            else:
                print(bcolors.GREY + "The move you entered consisted of more than " + str(maxNumMoves) + " cells"
                      + bcolors.ENDC)

        newCells.pop(0)  # remove the first cell, in order to only return the newly selected cells
        return newCells


# Main game loop function
def playGame(nimBoard, emptyCells, computerStrategy, firstTurn, playerLetter, computerLetter, maxNumMoves, N):
    # firstTurn indicates who shall play the first move (its possible values are 'computer' and 'player')
    playerTurn = False  # indicated if it is the HUMAN player's turn (if TRUE) or not (if FALSE)
    if firstTurn == "player":
        playerTurn = True

    playerLastMove = []  # contains the last move of the HUMAN player

    while not isBoardFull(nimBoard, N):  # run the game loop while there are still empty cells left on the 2D NIM board
        drawNimPalette(nimBoard, N)  # print current board state prior to any player move

        if playerTurn:
            cells = []
            playerLastMove.clear()  # clear the move of the previous turn, to store the move of the current turn
            # get and check the first cell of the HUMAN player's move
            firstCell = getAndCheckInputCells(nimBoard, "first", cells, maxNumMoves, N)
            cells.append(firstCell)  # append the first chosen cell to the cells to be player during the current turn
            writeCells(nimBoard, emptyCells, cells, playerLetter)  # "color" the cell using the player's letter (R/G)
            playerLastMove.append(firstCell)  # append the cell into the "last HUMAN player move" list

            # if the first cell is not a diagonal one, the player may enter up to two more
            if not isCellDiagonal(firstCell, N):
                newCells = getAndCheckInputCells(nimBoard, "next", cells, maxNumMoves,
                                                 N)  # contains the newly selected cells
                writeCells(nimBoard, emptyCells, newCells, playerLetter)
                playerLastMove = playerLastMove + newCells

        else:  # computer turn
            if len(emptyCells) > 1:
                print("Computer response :")
            playUsingStrategy = True
            # There are at most 5 empty cells -> If the board corresponds to a winning state, the computer must be able
            # to find the correct winning move.
            # Else if the current state corresponds to a losing state, continue playing user the user-defined strategy
            if (N * N - nimBoard[0] <= 5):  # below 5 empty cells, the computer knows how to play.
                winPossible, winningMove = getComputerWinningMove(emptyCells, N)

                if winPossible:  # COMPUTER can win -> winningMove contains the optimal move
                    writeCells(nimBoard, emptyCells, winningMove, computerLetter)
                    playUsingStrategy = False

            # The computer will play according to the user-selected strategy if there are move than 5 empty cells
            # OR there are at most 5 empty cells whose board arrangement corresponds to a losing state.
            if playUsingStrategy:
                match computerStrategy:
                    case "random":
                        chosenCells = getComputerMove_random(nimBoard, emptyCells, N)
                        writeCells(nimBoard, emptyCells, chosenCells, computerLetter)

                    case "first free":
                        chosenCells = getComputerMove_firstfit(nimBoard, emptyCells, N)
                        writeCells(nimBoard, emptyCells, chosenCells, computerLetter)

                    case "copycat":
                        chosenCells = getComputerMove_copycat(nimBoard, emptyCells, playerLastMove, N)
                        writeCells(nimBoard, emptyCells, chosenCells, computerLetter)

        # Complement the last active player flag, effectively causing an active player switch
        playerTurn = not playerTurn

    # Game finished -> If playerTurn = False then the last player who made a move was the HUMAN player, since the flag
    # got switched from True to False, prior to exiting the main game loop.
    # Else if playerTurn = True, that means that during the last turn playerTurn was set to False, i.e. the COMPUTER
    # was the last one to move.
    lastPlayer = 'Player' if not playerTurn == True else 'Computer'
    print(lastPlayer.upper() + " wins !!")


######### MAIN PROGRAM BEGINS #########
screen_clear()
print(bcolors.HEADER + """
---------------------------------------------------------------------
                     2-Dimensional NIM Game: RULES (I)
---------------------------------------------------------------------
    1.      A human PLAYER plays against the COMPUTER.
    2.      The starting position is an empty NxN board.
    3.      One player (the green) writes G, the other player 
                (the red) writes R, in empty cells.
""" + bcolors.ENDC)

input("Press ENTER to continue...")
screen_clear()

print(bcolors.HEADER + """
---------------------------------------------------------------------
                     2-Dimensional NIM Game: RULES (II) 
---------------------------------------------------------------------
    4.      The cells within the NxN board are indicated as 
            consecutive numbers, from 1 to N^2, starting from the 
            upper-left cell. E.g. for N=4, the starting position 
            and some intermediate position of the game would be 
            like those:
                    INITIAL POSITION        INTERMEDIATE POSITION
                    =====================   =====================
                    [  1 |  2 |  3 |  4 ]   [  1 |  2 |  3 |  4 ]
                    ---------------------   ---------------------
                    [  5 |  6 |  7 |  8 ]   [  5 |  R |  7 |  8 ]    
                    ---------------------   ---------------------
                    [  9 | 10 | 11 | 12 ]   [  9 |  R | 11 | 12 ] 
                    ---------------------   ---------------------
                    [ 13 | 14 | 15 | 16 ]   [  G |  G | 15 |  G ] 
                    =====================   =====================
                       COUNTER = [ 0 ]         COUNTER = [ 5 ]
                    =====================   =====================
""" + bcolors.ENDC)

input("Press ENTER to continue...")
screen_clear()

print(bcolors.HEADER + """
---------------------------------------------------------------------
                     2-Dimensional NIM Game: RULES (III) 
---------------------------------------------------------------------
    5.      In each round the current player's turn is to fill with 
            his/her own letter (G or R) at least one 1 and at most 
            3 CONSECUTIVE, currently empty cells of the board, all 
            of them lying in the SAME ROW, or in the SAME COLUMN 
            of the board. Alternatively, the player may choose ONLY
            ONE empty diagonal cell to play.
    6.      The player who fills the last cell in the board WINS.
    7.      ENJOY!!!
---------------------------------------------------------------------
""" + bcolors.ENDC)

maxNumMoves = 3

playNewGameFlag = True

while playNewGameFlag:
    if not startNewGame():
        break

    N = getBoardSize()
    emptyCells = [*range(1, (N * N) + 1)]  # list that contains the empty cells of the board
    nimBoard = initializeBoard(N)

    playerLetter, computerLetter = inputPlayerLetter()
    turn = whoGoesFirst()
    computerStrategy = howComputerPlays()

    print(bcolors.MSG + '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' + bcolors.ENDC)
    print(bcolors.MSG + 'A new ' + str(N) + 'x' + str(
        N) + ' game is about to start. The ' + turn + ' makes the first move.' + bcolors.ENDC)
    print(
        bcolors.MSG + ' * The computer will play according to the ' + bcolors.HEADER + computerStrategy + bcolors.MSG + ' strategy.' + bcolors.ENDC)
    print(
        bcolors.MSG + ' * The player will use the letter ' + playerLetter + ' and the computer will use the ' + computerLetter + '.' + bcolors.ENDC)
    print(bcolors.MSG + ' * The first move will be done by the ' + bcolors.BLUE + turn + '.' + bcolors.ENDC)
    print(bcolors.MSG + '---------------------------------------------------------------------' + bcolors.ENDC)

    # call the game playing function
    playGame(nimBoard, emptyCells, computerStrategy, turn, playerLetter, computerLetter, maxNumMoves, N)

    print(bcolors.MSG + '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' + bcolors.ENDC)
    print("Final board state : ")
    drawNimPalette(nimBoard, N)
######### MAIN PROGRAM ENDS #########
