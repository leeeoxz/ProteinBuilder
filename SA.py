#Leonardo Santos#
#June, 2017#

# Simulated Annealing #
import random
import copy
import math

class SimAnne:
	temperature = None #Initial temperature
	solution = None #Initial Solution
	dimension = None #Solution dimension
	lb = None #Lower boundery for solutions
	ub = None #Upper boundery for solutions
	lb_flag = None #Tells if exists a lower boundery for solutions
	ub_flag = None #Tells if exists a upper boundery for solutions
	temperature_decrease_factor = 0.9 #Decrease factor that will multiply the temperature 

	def __init__(self, dimension, lb=[], ub=[],temperature=1000):
		self.temperature = temperature
		self.lb = lb
		self.ub = ub
		self.dimension = dimension
		self.solution = self.create_solution() #creates the first solution with random numbers


	def create_solution(self):
		sol = []

		for d in xrange(self.dimension):
			sol.append(str(random.uniform(self.lb[d], self.ub[d]))) #Creates a solution with dimension components, ranging between bounderies
		return sol

	def getProbability (self, Ene1, Ene2, temperature):
		prob = math.exp(-((Ene2-Ene1)/self.temperature)) #Calculate the probability of some worst solution to be accepted
		return prob

	def updateTemperature (self):
		self.temperature = self.temperature*self.temperature_decrease_factor #Descrease the temperature by 0.9

	def getTemperature (self):
		return self.temperature
