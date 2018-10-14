from sympy import *
from itertools import combinations

DEBUG_OUTPUT = True

def debug_output(string):
    if DEBUG_OUTPUT:
        print(string)

def read_tasks(path : "path to the file with task def"):
    """Reads tasks in proprietary undocumented format,
       which is easily reversible, returns list of tasks"""
    lst = []
    with open(path, 'r') as file:
        for line in file:
            tmp = eval(line)
            lst.append(tmp)
    return(lst)        

def get_equation(task):
    substitutions = [(sympify(item[0]), sympify(item[1])) for item in task[2]]
    
    left  = sympify(task[0]).subs('C', 1)
    right = sympify(task[1]).subs('C', 1)
    
    tmp = right.subs([(sympify(item[0]), E) for item in task[2]])
    powers = ''.join(str(expand_log(ln(tmp), force = True)).split('+'))
    powers = ' '.join(powers.split())
    
    eq = (left.subs(substitutions),
          right.subs(substitutions),
          powers,
          left, 
          right)
    return eq

def linearize(equation, var_subs):
    tmp = (expand_log(ln(equation[0].subs(var_subs).simplify()).simplify(), 
                      force = True),
           expand_log(ln(equation[1].subs(var_subs).simplify()).simplify(), 
                      force = True))
    
    return tmp[0] - tmp[1]
    
def make_system(equation):
    lst = []
    lst.append(linearize(equation, [('L', 1), ('M', 1), ('T', E)]))
    lst.append(linearize(equation, [('L', 1), ('M', E), ('T', 1)]))
    lst.append(linearize(equation, [('L', E), ('M', 1), ('T', 1)]))
    return lst
    
def extract_vars(equation):
    """Separates free variables"""
    all_vars = equation[2].split()
    
    for subset in combinations(all_vars, 3):
        yield (symbols(list(set(all_vars) - set(subset))),
               symbols(list(subset)))
    
def optimagick(equation, solution, variables):
    """Does the main optimagick thing using good old logarithms
       and differentiation, fun stuff"""
    sol = None
    for item in solution:
         sol = item

    substitutions = [(variables[0][i], sol[i]) for i in range(0, len(sol))]
    tmp = expand_log(ln(equation[4].subs(substitutions)), force = True)
    
    substitutions = [(variables[1][i], 0) for i in range(0, len(sol))]
    yield (equation[3], exp(tmp.subs(substitutions)))
    
    for var in variables[1]:
        yield (equation[3], exp(diff(tmp, var)))

def solve_task(task):
    """Solves a single task, TBD: make latex"""
    equation = get_equation(task)
    debug_output("Equation after substitutions:")
    debug_output(equation[1])
    
    system = make_system(equation)
    debug_output("Linear system:")
    debug_output(system)
    
    for variables in extract_vars(equation):
        print("Trying to solve for these vars:" + str(variables))
        
        solution = linsolve(system, variables[0])
        debug_output(solution)
        
        if solution == EmptySet():
            debug_output("No solutions for these vars!")
            continue
        
        for eq in optimagick(equation, solution, variables):
            print(eq)

for task in read_tasks('task.txt'):
    solve_task(task)