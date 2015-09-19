#!/usr/bin/env python

__author__ = 'Batylan Nurbekov & Ari Goodman & Doruk Uzunoglu & Miguel Mora'

import sys, re, math, random, logging, time

DEBUG = 1
PUZZLE1_INIT_POP = 5
MU = .999

def getArgs():
    validArgNum = len(sys.argv) == 4

    if not validArgNum:
        raise Exception("The number of arguments should be 3")

    return (int(sys.argv[1]), sys.argv[2], int(sys.argv[3]))

class Puzzle(object):
    def __init__(self, filePath, secondsToWork):
        self.filePath = filePath
        self.secondsToWork = secondsToWork
        self.population = set()
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

    #Creates a gene with the specified length with random values in range [min, max]
    def createGene(self, min, max, length, factorial):
        if not factorial:
            return tuple(random.randint(min, max) for x in range(length))
        else:
            return tuple(random.randint(min, max) for x in range(length)) #TODO


    def weighted_choice(self, weights):
        total = sum(weights)
        r = random.uniform(0, total)
        upto = 0
        index = 0
        for w in weights:
            if upto + w >= r:
                return index

            upto += w
            index += 1

    def competition(self, a, b):
        return

    #Abstract methods:
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

        #if the number you are paying with is bigger than the number of possible states, use the number of possible states
        if PUZZLE1_INIT_POP > math.pow(2, len(self.input)):
            initPopSize = math.pow(2, len(self.input))
        else:
            initPopSize = PUZZLE1_INIT_POP


        for i in range(initPopSize):
            self.population.add(self.createGene(0, 1, len(self.input), True))

        logging.debug("Init population %s", self.population)

        self.findSolution()

        logging.debug("Solution %s", self.solution)

    def findSolution(self):
        while True:
            # Check timer and convergence TODO: add convergence check
            if self.timerElapsed():
                break

            # Filter out the duplicates (converges to solution slower, but creates diversity in population)
            # TODO: check the implementation with duplicates (the list instead of set will have to be implemented for pop)

            self.fitness = []
            bestFitness = 0

            # Estimate fitness for each gene & return solution, if found
            for gene in self.population:
                fitness = self.estimateFitness(gene)

                if fitness == self.goal:
                    self.solution = gene
                    return

                if fitness > bestFitness:
                    bestFitness = fitness
                    self.solution = gene

                self.fitness.append(fitness)

            logging.debug("Fitness: %s", self.fitness)

            # self.winners = []
            # # Pick genes based on fitness
            # for i in range(0,len(self.fitness),2):
            #     self.winners.append(self.competition(self.fitness[i],self.fitness[i+1]))


            # Perform crossover

            # Perform mutation

    def estimateFitness(self, gene):
        sum =0
        i = 0
        for chromosome in gene:
            sum+= chromosome*self.input[i]
            i+=1
        difference = sum-self.goal
        if difference <= 0:
            return sum
        else:
            return 1/difference

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
