# EuclidGCD
Python3 module including functions to compute GCD of integer numbers (i.e. set Z), even modular, or polynomials over a field (e.g. R, Z3[x])

The module also includes helper functions to compute division or multiplication between polynomials.

It can also be used as a standalone program from the terminal.

Usage:

EuclidMCD.py [-h] [-d] [--mod-inverse] [--poly] [--poly-div] [--poly-mul] vals vals [modulus]

Euclid's algorithm to compute GCD

positional arguments:

  - vals -> integers or comma-separated list of integers (representing polynomial's coefficients with decreasing degree order). 
  
    **Examples:**
    - **python3 EuclidMCD.py 6 15** => (3, 1, -2).  
    The GCD is the first element of the tuple. The remaining 2 elements are the coefficients of a linear combination of the given numbers which results equal to their GCD.  
    In this case, GCD(6, 15) = 3 = 1 * 15 - 2 * 6.  
    ***NOTE:** the second element of the tuple is **ALWAYS** the coefficient for max(a, b), where a and b are the values passed to the program (6 and 15 in the example). Obviously, this means that the third element is **ALWAYS** the coefficient for min(a,b).*
    
    - **python3 EuclidMCD.py 1,0,-1 1,1 --poly** => 1.0 x + 1.0.  
    The values passed to the program are interpreted as the coefficients of 2 polynomials, so it computes GCD(x^2 - 1; x + 1) = x + 1. 
- modulus (optional) -> integer

optional arguments:

  - -h, --help  
  show an help message and exit
  
  - -d, --debug  
  Enables debug printing, which prints all the passages of the Euclid's algorithm applied to arguments
  
  - --mod-inverse  
  Computes the modular inverse of the first positional argument using the second positional argument as modulus.  
  (e.g. **python3 EuclidMCD.py 3 5 --mod-inverse** => 2)
  
  - --poly  
  Interpret the positional arguments as coefficients of 2 polynomials in the real field, and compute the GCD between them
  
  - --poly-div  
  Interpret the positional arguments as coefficients of 2 polynomials in the real field, and compute the division between them 
  
  - --poly-mul  
  Interpret the positional arguments as coefficients of 2 polynomials in the real field, and compute the multiplication between them
  
  ***NOTE:** options --poly, --poly-div and --poly-mul can be combined with the positional argument **modulus** in order to perform modular computations*
