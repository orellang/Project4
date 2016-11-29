# Program to solve Travelling Salesman Problem
# Uses portions of tsp-verifier.py provided by instructor
#!/usr/bin/python

import math
import sys
import re	

# Loads list of cities and coordinates from file and stores in dictionary
# Takes file name as parameter
# Returns dictionary with cities and coordinates
def loadInput(fileName):

        # Creates a dictionary of cities from file with x and y coord as values
        citiesList = {}
        with open(fileName) as file:
                for line in file:
                        lineparse = re.findall(r'[^,;\s]+', line)
                        citiesList[int(lineparse[0])] = (int(lineparse[1]),int(lineparse[2]))
                        #citiesList.append([int(lineparse[0]),int(lineparse[1]),int(lineparse[2])])
                
        return citiesList

# Calculates euclidean distance between two cities and rounds to nearest integer
# Takes two cities with x and y locations as index 0, 1 respectively
# Returns distance between cities
def getDistance(a, b):
        # a and b are sets of three integers with their id, x and y coordinates as
        # Calculates the Euclidean distance between points, rounded to nearest int
        
        dx = a[0]-b[0]
        dy = a[1]-b[1]

        #Calculate distance using Pythagorean Theorems and round to nearest int
        return int(round(math.sqrt(dx*dx + dy*dy)))

# Creates matrix with distances from each city to every other city
# Takes dictionary with cities and coordinates as parameter
# Returns distance matrix
def createDistanceMatrix(cities):
        # Creates a matrix of the distance between all cities
        matrix = [[0]*len(cities) for _ in range(len(cities))]
        for i in range(len(cities)):
                for j in range(i+1, len(cities)):
                        matrix[i][j] = matrix[j][i] = getDistance(cities[i], cities[j])
        return matrix   

# Finds minimum spanning tree for cities using Prim's algorithm
# Handles cases with cities in same location
# Takes dictionary of cities, distance matrix, and starting city as parameters
# Returns minimum spanning tree
def MST(cities, matrix, start): 

        vertices = {}
        mst = {}
        for city in cities:
                vertices[city] = [start, matrix[start][city]]
                mst[city] = []
        #print(mst)
        #print(len(cities))
        i = 0
        size = len(cities) - 1

        while i < size:
                shortestDist = float("inf")
                nearestCity = None
                for v in vertices:
                        if vertices[v][1] > 0 and vertices[v][1] < shortestDist:
                                shortestDist = vertices[v][1]
                                nearestCity = v
                for v in vertices:
                        if cities[v][0] == cities[nearestCity][0] and cities[v][1] == cities[nearestCity][1]:
                                if v != nearestCity:
                                        vertices[v][1] = 0
                                        mst[v].append(vertices[v][0])
                                        mst[vertices[v][0]].append(v)
                                        i += 1
                vertices[nearestCity][1] = 0
                mst[nearestCity].append(vertices[nearestCity][0])
                mst[vertices[nearestCity][0]].append(nearestCity)
                del vertices[nearestCity]
                
                for v in vertices:
                        if vertices[v][1] > matrix[v][nearestCity]:
                                vertices[v][1] = matrix[v][nearestCity]
                                vertices[v][0] = nearestCity
                i += 1

        return mst

# Performs depth first search on graph
# Takes graph and starting vertex as parameters
# Returns adjacency list in dfs order
def dfs(graph, start):
        visited = []
        stack = [start, ]

        while stack:
                node = stack.pop()
                
                if node not in visited:
                        visited.append(node)
                        stack.extend([x for x in graph[node] if x not in visited])
        return visited

# Finds total distance of tour between cities
# Takes list of cities in order visited and dictionary of cities with coordinates as parameters
# Returns distance traveled
def tourDistance(tour, cities):
        n = len(cities)
        dist = 0
        for i in range(n):
                dist = dist + getDistance(cities[tour[i]],cities[tour[i - 1]])
        return dist

        
# Finds vertices with odd degrees in minimum spanning tree
# Takes mst as parameter
# Returns list of odd vertices
def oddDegVert(mst):
        oddVerts = []
        for i in range (len(mst)):
                
                if len(mst[i])%2 == 1:
                        oddVerts.append(i)
        return oddVerts

# Finds minimal matching for odd degree vertices and adds edges
# Takes list of odd degree vertices, mst, and distance matrix as parameters
# Modifies mst parameter by adding edges
def minimalMatching(oddVerts, mst, matrix):
        closestVert = None
        while oddVerts:
                current = oddVerts.pop()
                shortestDist = float("inf")
                for v in oddVerts:
                        if matrix[current][v] < shortestDist:
                                shortestDist = matrix[current][v]
                                closestVert = v
                mst[current].append(closestVert)
                mst[closestVert].append(current)
                oddVerts.remove(closestVert)

# Calculates eulerian tour of graph
# takes graph and starting vertex as parameters
# returns eulerian tour path as adjacency list of vertices
def eulerianTour(graph, start):
        path = []
        stack = []
        stack.append(start)

        while stack:
                current = stack[len(stack)-1]
                if len(graph[current]) > 0:
                        stack.append(graph[current][0])
                        graph[graph[current][0]].remove(current)
                        del graph[current][0]
                else:
                        path.append(stack.pop())
        return path

# finds hamiltonian path from eulerian tour
# takes eulerian tour in form of adjacency list as parameter
# returns hamiltonian path as list of vertices
def findHamPath(eulerianTour):
        
        visited = {}
        hamPath = []
        for vertex in eulerianTour:
                visited[vertex] = 0
                
        start = eulerianTour[0]
        visited[start] = 1
        hamPath.append(start)

        for vertex in eulerianTour:
                if visited[vertex] == 0:
                        hamPath.append(vertex)
                        visited[vertex] = 1
        return hamPath

# Writes tsp approximation solution to file
# takes original problem filename, distance of tsp path, and the path as parameters
# writes the solution to file named the original filename with ".tour" appended
# Solution file is formatted with distance on first line and city IDs on each following line
def writeSolution(filename, distance, path):
        outputFile = filename + ".tour"
        output = open(outputFile, 'w+')
        output.write(str(distance)+"\n")
        for i in path:
                output.write(str(i)+"\n")        

# Optimizes hamiltonian path using 2-Opt method
# Takes hamiltonian path, dictionary of cities and coordinates, and the number of times to
# attempt optimization as paramters
# Calls twoOpt function until number of tries has been reached or until solution is no longer improved. 
# Returns optimized path as list of vertices
def optimize(hamPath, cities, tries):
        changeMade = True
        optimizedPath = list(hamPath)
        optimizedDistance = tourDistance(optimizedPath, cities)
        i = 0
        while changeMade == True and i < tries:
                changeMade = False
                newPath = twoOpt(optimizedPath, cities)
                newDistance = tourDistance(newPath, cities)
                if newDistance < optimizedDistance:
                        optimizedPath = newPath
                        changeMade = True
                        i += 1
        #print(optimizedPath)
        return optimizedPath

# Takes hamiltonian path and tries swapping edges to see if path is improved
# Takes path and dictionary of cities with coordinates as parameters
# returns best path found as list of vertices
def twoOpt(hamPath, cities):

        changeMade = True
        bestPath = list(hamPath)
        bestDistance = tourDistance(hamPath, cities)
        for i in range(0, len(bestPath)-2):
                for j in range(i+1, len(bestPath)-1):
                        newPath = twoOptSwap(bestPath, i, j)
                        newDistance = tourDistance(newPath, cities)
                        if newDistance < bestDistance:
                                bestPath = newPath
        return bestPath

# Swaps edges in path given two vertices
# Takes path and vertices to swap as parameters
# Swaps order of vertices between two chosen vertices
# Returns new path as list of vertices
def twoOptSwap(path, i, j):
        newPath = list(path)
        newPath[i:j] = newPath[i:j][::-1]
        return newPath

# Finds tsp approximation from problem file and writes solution to file
# Takes problem filename as parameter
# Writes solution to file and prints tour and distance to console
def tsp(filename):

        cities = loadInput(filename)
        matrix = createDistanceMatrix(cities)
        tour = MST(cities, matrix, 0)
        shortTour = dfs(tour, 0)
        oddVerts = oddDegVert(tour)
        minimalMatching(oddVerts, tour, matrix)
        euTour = (eulerianTour(tour, 0))
        #print(tourDistance(euTour,cities))
        hamPath = findHamPath(euTour)
        optPath = optimize(hamPath, cities, 1)
        #print("Opt Path: ")
        print(optPath)
        distance = tourDistance(optPath, cities)
        print(distance)
        writeSolution(filename, distance, optPath)
        
                

sys.argv = ["tsp.py", "tsp_example_1.txt"]
def main(argv):
        filename = sys.argv[1]
        tsp(filename)


if __name__ == "__main__":
        main(sys.argv[1])
