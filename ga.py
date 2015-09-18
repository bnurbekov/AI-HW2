#!/usr/bin/env python

__author__ = 'Batylan Nurbekov & Ari Goodman & Doruk Uzunoglu & Miguel Mora'

import sys, re, math, random, logging, time

DEBUG = 1
PUZZLE1_INIT_POP = 5

def getArgs():
    validArgNum = len(sys.argv) == 4

    if not validArgNum:
        raise Exception("The number of arguments should be 3")

    return (int(sys.argv[1]), sys.argv[2], int(sys.argv[3]))

class Puzzle(object):
    def __init__(self, filePath, secondsToWork):
        self.filePath = filePath
        self.secondsToWork = secondsToWork
        self.population = []
        self.input = []
        self.goal = None
        self.executionStart = None
        self.solution = None

    #Checks if the timer has elapsed
    def timerElapsed(self):
        elapsed = False

        if self.executionStart == None:
            self.executionStart = time.time()
        else:
            if (time.time() - self.executionStart) > self.secondsToWork:
                elapsed = True

        return elapsed

    def run(self):
        raise Exception("run() is not implemented.")

    def findSolution(self):
        raise Exception("findSolution() is not implemented.")

    def estimateFitness(self, gene):
        raise Exception("estimateFitness() is not implemented.")

    def parseFile(self):
        raise Exception("parseFile() is not implemented.")

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

class PuzzleOne(Puzzle):
    def run(self):
        logging.debug("Running puzzle one...")

        self.parseFile()

        if PUZZLE1_INIT_POP > math.pow(2, len(self.input)):
            initPopSize = math.pow(2, len(self.input))
        else:
            initPopSize = PUZZLE1_INIT_POP


        samplingList = [0, 1] * len(self.input)

        for num in self.input:
            self.population.append(random.sample(samplingList, initPopSize))

        logging.debug("Init population %s", self.population)

        self.findSolution()

        logging.debug("Solution %s", self.solution)

    def findSolution(self):
        while True:
            # Check timer
            if self.timerElapsed():
                break

            # Filter out the duplicates TODO: Ask Beck about duplicates

            # Estimate fitness for each gene & return solution, if found

            # Pick genes based on fitness

            # Perform crossover

            # Perform mutation


    def estimateFitness(self, gene):
        raise Exception("estimateFitness() is not implemented.")

    def parseFile(self):
        file = open(self.filePath, "r")
        rows = file.readlines()

        self.input = []
        count = 0

        for row in rows:
            chars = filter(None, re.split('\t|\s|\n|\v|\r', row))
            num = int(chars[0])

            if count == 0:
                self.goal = num
            else:
                self.input.append(num)

            count += 1

class PuzzleTwo(Puzzle):
    def run(self):
        return

    def parseFile(self):
        return

class PuzzleThree(Puzzle):
    def run(self):
        return

    def parseFile(self):
        return

if __name__ == "__main__":
    (puzzleNum, filePath, secs) = getArgs()

    logging.basicConfig(level=logging.DEBUG)

    puzzle = PuzzleFactory.initPuzzle(puzzleNum, filePath, secs)
    puzzle.run()
