__author__ = 'Batylan Nurbekov & Ari Goodman & Doruk Uzunoglu & Miguel Mora'

import sys

def getArgs():
    validArgNum = len(sys.argv) == 4

    if not validArgNum:
        raise Exception("The number of arguments should be 3")

    return (int(sys.argv[1]), sys.argv[2], int(sys.argv[3]))

class PuzzleFactory:
    @staticmethod
    def initPuzzle(puzzleNum, filePath, secondsToWork):
        if (puzzleNum == 1):
            return PuzzleOne(filePath, secondsToWork)
        elif (puzzleNum == 2):
            return PuzzleTwo(filePath, secondsToWork)
        elif (puzzleNum == 3):
            return PuzzleThree(filePath, secondsToWork)
        else:
            return None

class Puzzle(object):
    def __init__(self, filePath, secondsToWork):
        self.filePath = filePath
        self.secondsToWork = secondsToWork

    def run(self):
        raise Exception("Run() is not implemented.")

class PuzzleOne(Puzzle):
    def run(self):
        return

class PuzzleTwo(Puzzle):
    def run(self):

        return

class PuzzleThree(Puzzle):
    def run(self):

        return

if __name__ == "__main__":
    (puzzleNum, filePath, secs) = getArgs()

    puzzle = PuzzleFactory.initPuzzle(puzzleNum, filePath, secs)
    puzzle.run()
