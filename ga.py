#!/usr/bin/env python

__author__ = 'Batylan Nurbekov & Ari Goodman & Doruk Uzunoglu & Miguel Mora'

import sys, re, math, random, logging, time

DEBUG = 1
PUZZLE1_INIT_POP = 50
PUZZLE2_INIT_POP = 6
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
        self.population = {}
        self.input = []
        self.goal = None
        self.executionStart = None
        self.solution = None
        self.bestFitness = -1

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


    def weighted_choice(self, population):
        total = sum(weight for gene, weight in population.iteritems())
        r = random.uniform(0, total)
        upto = 0

        for gene, weight in population.iteritems():
            if upto + weight >= r:
                return gene

            upto += weight

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
            self.population[self.createGene(0, 1, len(self.input), False)] = 0

        self.initPopSize = len(self.population)

        logging.debug("Init population %s", self.population)

        self.findSolution()

        logging.debug("Solution %s", self.solution)

    def findSolution(self):
        while True:
            # Check timer and convergence TODO: add convergence check
            if self.timerElapsed():
                break

            # Estimate fitness for each gene & return solution, if found
            for gene in self.population:
                fitness = self.estimateFitness(gene)

                #Solution found, exit
                if fitness == self.goal:
                    self.solution = gene
                    return

                #Tries to estimate the best fitness so far
                if fitness > self.bestFitness:
                    self.bestFitness = fitness
                    self.solution = gene

                self.population[gene] = fitness

            logging.debug("Population: %s", self.population)

            # Perform crossover
            children = {}

            while True:
                if self.initPopSize < len(children) + 2:
                    break

                parent1 = self.weighted_choice(self.population)
                parent2 = self.weighted_choice(self.population)

                #generate cut-off (split)
                split = random.randint(1, len(parent1)-1)
                child1 = parent1[0:split] + parent2[split:len(parent1)]
                child2 = parent2[0:split] + parent1[split:len(parent1)]

                # Perform mutation
                if random.uniform(0, 1) < 0.05:
                    child1lst = list(child1)
                    child1lst[random.randint(0, len(child1lst)-1)] ^= 1
                    child1 = tuple(child1lst)

                if random.uniform(0, 1) < 0.05:
                    child2lst = list(child2)
                    child2lst[random.randint(0, len(child2lst)-1)] ^= 1
                    child2 = tuple(child2lst)

                children[child1] = 0
                children[child2] = 0

            self.population = children

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
        logging.debug("Running puzzle two...")

        self.parseFile()

        #if the number you are paying with is bigger than the number of possible states, use the number of possible states
        if PUZZLE2_INIT_POP > math.pow(2, len(self.input)):
            initPopSize = math.pow(2, len(self.input))
        else:
            initPopSize = PUZZLE2_INIT_POP

        for i in range(initPopSize):
            self.population.add(self.createGene(0, 2, len(self.input), 2))

        logging.debug("Init population %s", self.population)

        self.findSolution()

        logging.debug("Solution %s", self.solution)
        return

    def parseFile(self):
        file = open(self.filePath, "r")
        rows = file.readlines()

        self.input = []
        count = 0

        for row in rows:
            chars = filter(None, re.split('\t|\s|\n|\v|\r', row))
            num = int(chars[0])
	    self.input.append(num)
            count += 1

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
