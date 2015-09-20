#!/usr/bin/env python

__author__ = 'Batylan Nurbekov & Ari Goodman & Doruk Uzunoglu & Miguel Mora'

import sys, re, math, random, logging, time

DEBUG = 1
PUZZLE1_INIT_POP = 10
PUZZLE2_INIT_POP = 6
PUZZLE3_INIT_POP = 3
MU = .999

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
        self.population = {}
        self.input = []
        self.goal = None
        self.executionStart = None
        self.solution = None
        self.bestFitness = -1

        #Stat variables
        self.generationNum = 0
        self.generationWhenSolutionFound = 0

    #Runs puzzle
    def run(self):
        self.parseFile()

        initPopSize = self.getInitPopSize()
        maxPopSize = self.getMaxPopSize()

        #if the number you are paying with is bigger than the number of possible states, use the number of possible states
        if initPopSize > maxPopSize:
            initPopSize = maxPopSize

        for i in range(initPopSize):
            self.population[self.createGene()] = 0

        self.initPopSize = len(self.population)

        logging.debug("Init population %s", self.population)

        self.findSolution()

        logging.debug("Solution %s", self.solution)

        self.printStats()

    def findSolution(self):
        while True:
            # Check timer and convergence (including if solution is found)
            if self.timerElapsed() or self.converged():
                break

            #Update stats
            self.generationNum += 1

            # Estimate fitness for each gene
            for gene in self.population:
                fitness = self.estimateFitness(gene)

                #Tries to estimate the best fitness so far
                if fitness > self.bestFitness:
                    self.solution = gene
                    self.generationWhenSolutionFound = self.generationNum
                    self.bestFitness = fitness

                self.population[gene] = fitness

            logging.debug("Population: %s", self.population)

            children = {}

            while True:
                if self.initPopSize < len(children) + 2:
                    break

                parent1 = self.weighted_choice(self.population)
                parent2 = self.weighted_choice(self.population)

                #generate cut-off (split)
                split = random.randint(1, len(parent1)-1)

                # Perform crossover
                child1_lst = list(parent1[0:split] + parent2[split:len(parent1)])
                child2_lst = list(parent2[0:split] + parent1[split:len(parent1)])

                # Repair invalid genes (only applicable to Puzzle 2)
                self.repair(child1_lst)
                self.repair(child2_lst)

                # Perform mutation
                self.mutate(child1_lst)
                self.mutate(child2_lst)

                children[tuple(child1_lst)] = 0
                children[tuple(child2_lst)] = 0

            self.population = children

    #Checks if the timer has elapsed
    def timerElapsed(self):
        elapsed = False

        if self.executionStart == None:
            self.executionStart = time.time()
        else:
            if (time.time() - self.executionStart) > self.secondsToWork:
                elapsed = True

        return elapsed

    #Picks a sample from the population with
    def weighted_choice(self, population):
        total = sum(weight for gene, weight in population.iteritems())
        r = random.uniform(0, total)
        upto = 0

        for gene, weight in population.iteritems():
            if upto + weight >= r:
                return gene

            upto += weight

    #Parses file to extract input and goal representations
    def parseFile(self):
        file = open(self.filePath, "r")
        rows = file.readlines()

        data = []

        for row in rows:
            str_list = filter(None, re.split('\t|\s|\n|\v|\r', row))
            representation = self.getInputRepresentation(str_list)
            data.append(representation)

        self.processData(data)

    def printStats(self):
        print "Solution: "
        print self.convertRepresentationToString(self.solution)
        print "Solution score: %d" % self.getScore(self.solution)
        print "Generations the algorithm ran for: %d" % self.generationNum
        print "Generation when solution was found: %d" % self.generationWhenSolutionFound

    # Processes the input representation data extracted from file
    def processData(self, data):
        raise Exception("processData() is not implemented.")

    # Converts input string list into representation of an input used in the puzzle
    def getInputRepresentation(self, str_list):
        raise Exception("getInputRepresenation() is not implemented.")

    # Creates gene
    def createGene(self):
        raise Exception("createGene() is not implemented.")

    # Repairs a gene to conform the constraints of the puzzle
    def repair(self, gene_lst):
        raise Exception("repair() is not implemented.")

    # Mutates a gene
    def mutate(self, gene_lst):
        raise Exception("mutate() is not implemented.")

    # Estimates fitness of a gene
    def estimateFitness(self, gene):
        raise Exception("estimateFitness() is not implemented.")

    # Gets maximum population size
    def getMaxPopSize(self):
        return int(math.pow(2, len(self.input)))

    # Gets initial population size
    def getInitPopSize(self):
        raise Exception("getInitPopSize() is not implemented.")

    # Converts internal input representation defined by the gene to a string
    def convertRepresentationToString(self, gene):
        raise Exception("convertRepresentationToString() is not implemented.")

    def getScore(self, gene):
        raise Exception("getScore() is not implemented.")

    def converged(self):
        raise Exception("converged() is not implemented.")

class PuzzleOne(Puzzle):
    def getInputRepresentation(self, str_list):
        return float(str_list[0])

    def processData(self, data):
        self.goal = data[0]
        self.input = data[1:len(data)]

    def createGene(self):
        return tuple(random.randint(0, 1) for x in range(len(self.input)))

    def getInitPopSize(self):
        return PUZZLE1_INIT_POP

    def converged(self):
        return self.bestFitness == self.goal

    def repair(self, gene_lst):
        #No need to repair genes for this puzzle
        return

    def mutate(self, gene_lst):
        if random.uniform(0, 1) < 0.05:
            gene_lst[random.randint(0, len(gene_lst)-1)] ^= 1

    def estimateFitness(self, gene):
        sum = 0
        i = 0

        for chromosome in gene:
            sum += chromosome * self.input[i]
            i += 1

        difference = sum - self.goal

        if difference <= 0:
            return sum
        else:
            return 1/difference

    def convertRepresentationToString(self, gene):
        return ' '.join([str(self.input[i]) for i in range(len(gene)) if gene[i]])

    def getScore(self, gene):
        score = sum([self.input[i] for i in range(len(gene)) if gene[i]])

        return score if score <= self.goal else 0

class PuzzleTwo(Puzzle):
    def getInputRepresentation(self, str_list):
        return float(str_list[0])

    def processData(self, data):
        self.input = data

    def createGene(self):
        available_bins = {0: 0, 1: 0, 2: 0}
        gene = []

        for i in range(len(self.input)):
            bin_index = random.choice(available_bins.keys())
            gene.append(bin_index)
            available_bins[bin_index] += 1
            if available_bins[bin_index] >= len(self.input)/3:
                del available_bins[bin_index]

        return tuple(gene)

    def getInitPopSize(self):
        return PUZZLE2_INIT_POP

    def repair(self, gene_lst):
        bins = {0:[], 1:[], 2:[]}

        while True:
            i = 0
            for chromosome in gene_lst:
                bins[chromosome].append(i)
                i += 1

            if (len(bins[0]) == len(bins[1])) and (len(bins[1]) == len(bins[2])):
                break

            #find bins with largest number of elements and swap elements that are the same in both lists
            largest_bin = max(bins.iterkeys(), key=(lambda key: len(bins[key])))
            smallest_bin = min(bins.iterkeys(), key=(lambda key: len(bins[key])))

            rand_index_in_largest_bin = random.choice(bins[largest_bin])

            gene_lst[rand_index_in_largest_bin] = smallest_bin

            for bin_i in bins.iterkeys():
                bins[bin_i] = []

    def mutate(self, gene_lst):
        if random.uniform(0, 1) < 0.05:
            available_bins = [0, 1, 2]

            #Select a bin randomly
            firstBin = random.choice(available_bins)
            available_bins.remove(firstBin)

            #Select a random number in that bin
            a = len(self.input)/3 - 1
            index_in_bin1 = random.randint(0, len(self.input)/3 - 1)
            count = 0
            for i in range(len(gene_lst)):
                if gene_lst[i] == firstBin:
                    if count == index_in_bin1:
                        index_in_gene1 = i
                        break

                    count += 1

            #Select another bin randomly
            secondBin = random.choice(available_bins)

            #Select a random number in that bin
            index_in_bin2 = random.randint(0, len(self.input)/3 - 1)
            count = 0
            for i in range(len(gene_lst)):
                if gene_lst[i] == secondBin:
                    if count == index_in_bin2:
                        index_in_gene2 = i
                        break

                    count += 1

            #Switch the two random numbers found
            temp = gene_lst[index_in_gene1]
            gene_lst[index_in_gene1] = gene_lst[index_in_gene2]
            gene_lst[index_in_gene2] = temp

    def estimateFitness(self, gene):
        score = self.getScore(gene)

        if (score) > 0:
            return score
        else:
            return 1/math.fabs(score)

    def convertRepresentationToString(self, gene):
        return '\n'.join([' '.join(["Bin #1:"]+[str(self.input[i]) for i in range(len(gene)) if gene[i] == 0]),
                         ' '.join(["Bin #2:"]+[str(self.input[i]) for i in range(len(gene)) if gene[i] == 1]),
                         ' '.join(["Bin #3:"]+[str(self.input[i]) for i in range(len(gene)) if gene[i] == 2])])

    def getScore(self, gene):
        product = 1
        sum = 0

        i = 0
        for chromosome in gene:
            if chromosome == 0:
                product *= self.input[i]
            elif chromosome == 1:
                sum += self.input[i]

            i += 1

        return (product + sum) / 2

    def converged(self):
        return False

class PuzzleThree(Puzzle):
    def getInputRepresentation(self, string):
        return

    def createGene(self):
        return

    def getInitPopSize(self):
        return PUZZLE3_INIT_POP

    def repair(self, gene_lst):
        #No need to repair genes for this puzzle
        return

    def mutate(self, gene_lst):
        return

    def estimateFitness(self, gene):
        return

    def convertRepresentationToString(self, gene):
        return

    def getScore(self, gene):
        return

    def converged(self):
        return False

if __name__ == "__main__":
    (puzzleNum, filePath, secs) = getArgs()

    logging.basicConfig(level=logging.DEBUG)

    puzzle = PuzzleFactory.initPuzzle(puzzleNum, filePath, secs)
    puzzle.run()
