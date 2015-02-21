import sys
from eigenproblem_notebook import Problem_Notebook, Eigenproblem

class Menu:
  '''Display a menu and respond to choices when run.'''
  def __init__(self):
    self.notebook = Problem_Notebook()
    self.choices = {
      1: self.show_problems,
      2: self.add_problem,
      3: self.quit
    }

  def display_menu(self):
    print("""
      Notebook Menu
      1. Show all Problems
      2. Add Problem
      3. Quit 
    """)

  def run(self):
    '''Display the menu and respond to choices.'''
    while True:
      self.display_menu()
      choice = input("Enter an option: ")
      action = self.choices.get(choice)
      if action:
        action()
      else:
        print("{0} is not a valid choice".format(choice))

  def show_problems(self, notes=None):
    if not notes:
      problems = self.Notebook.problems
      for problem in problems:
        print("{0}: {1}\n{2}".format(
          problem.id, problem.creation_date, problem.dimension))
  
  def quit(self):
    print("Thank you for using your notebook today.")
    sys.exit(0)

if __name__ == "__main__":
  m = Menu()
  m.run()