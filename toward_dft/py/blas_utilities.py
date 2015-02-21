from numpy import *
from numpy.random import rand
import datetime

class Eigenproblem:
  '''Represent an eigenproblem to be solved. 
  Typically we will be given an Hermitian operator and we will solve it by various methods.'''

  def __init__(self, dimension):
    '''Initialize an eigenproblem with a specific dimension.'''
    self.creation_date = datetime.date.today()
    self.dimension = dimension

  def random_operator(self):
    rando = rand(self.dimension, self.dimension)
    self.operator = rando

  def print_operator(self):
    print self.operator
    

class blas():
    def main(self):     
        return 

    def random_matrix(self, rows, cols):
      print("Build a random_matrix")
      A = rand(rows,cols)

      return A

    def print_matrix(rows, cols):
      print("Printing a matrix ")
      return

    def print_vector(n):
      return

    def random_vector(n):
      return

    def identity_matrix(n):
      return