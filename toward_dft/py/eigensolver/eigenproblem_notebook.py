from numpy import *
from numpy.random import rand
import datetime

last_id = 0

class Eigenproblem:
  '''Represent an eigenproblem to be solved. 
  Typically we will be given an Hermitian operator and we will solve it by various methods.'''

  def __init__(self, dimension):
    '''Initialize an eigenproblem with a specific dimension.'''
    self.creation_date = datetime.date.today()
    self.dimension = dimension
    self.vectors = []
    global last_id
    last_id += 1
    self.if = last_id

  def random_operator(self):
    rando = rand(self.dimension, self.dimension)
    self.operator = rando

  def random_vector(self):
    rando = rand(self.dimension)
    self.vectors.append(rando)

  def print_operator(self):
    print self.operator

  def print_vector(self, i=None):
    if i == None: print self.vectors
    else: print self.vectors[i]

  def identity_matrix(n):
      return

class Problem_Notebook(object):
  """Initialize a notebook with an empty list."""
  def __init__(self):
    self.problems = []

  def new_problem(self, dimension):
    self.problems.append(Eigenproblem(dimension))
  

if __name__ == "__main__":
  eig = Eigenproblem(4)
  eig.random_operator()
  eig.print_operator()
  eig.random_vector()
  eig.print_vector()