############################################################################################
#  Project Name:                    Bernoulli and Eular Symbols                            #
#  By                               Lin JIU                                                #
#  Project Initiate:                July 1st, 2018                                         #
#  Last Update:                     February 3rd, 2020                                     #
############################################################################################

############################################################################################
##    Section I: Evaluating polynomials by Bernoulli or Euler (umbral) symbols            ##
############################################################################################

def EvalSym(expression, variableList, degreeList, sym='B'):
    '''
    Given an 

    expression (mainly polynomial), 
    variablist = the list of variables,
    degreeList = the list of highest degrees of each variable,
    
    and

    sym, which can be one of the following options:
        'B', the default option    : evaluating by the Bernoulli symbol, namely, 
                B^k = B_k, the kth Bernoulli number;
        'E', the Euler symbol      : evaluating by the Euler symbol, namely, 
                E^k = E_k, the kth Euler number;
        'C', the Lucas' convertion : 
                C^k = B_k/k, the kth Bernoulli number divided by k. 

    For example, given a polynomial

    P(x, y, z) = 3 * x^2 * y - y^4 * z^2 + 4 * x * y * z -10, 

    the input, besides expression = P, should also include:

    1. variablelist = [x, y, z]
    2. degreeList = [2, 4, 2], since the highest degrees for x, y, and z are 
                2, 4, and 2, respectively. 

    Therefore,
    
    1. if one evaluate it by the Bernoulli symbol, 
          EvalSym(P, [x, y, z], [2, 4, 2]) 
        = EvalSym(P, [x, y, z], [2, 4, 2], 'B') 
        = 3 * B_2 * B_1 - B_4 * B_2 + 4 * (B_1)^3 - 10
        = 3 * (1 / 6) * (-1 / 2) - (-1 / 30) * (1 / 6) + 4 * (-1 / 2)^3 - 10 
        = -967 / 90;

    2. if evaluating by the Euler symbol, 
          EvalSym(P, [x, y, z], [2, 4, 2], 'E') 
        = 3 * E_2 * E_1 - E_4 * E_2 + 4 * (E_1)^3 - 10
        = 3 * (-1) * 0 - 5 * (-1) + 4 * 0^3 - 10 
        = -5.
    '''

    re = expression;                                                      # This is the re(turn). #
    n = len(variableList);                                                # the number of variables #
    for i in range(n):                                                    # Evaluate the variables, one by one #
        for j in range(degreeList[i])[::-1]:                              # The [::-1] part is to use the decreassing order #
            if sym == 'B':
                re = re.subs(variableList[i]^(j + 1) == bernoulli(j + 1))
            elif sym == 'E':
                re = re.subs(variableList[i]^(j + 1) == euler_number(j + 1))
            elif sym == 'C':
                re = re.subs(variableList[i]^(j+1) == bernoulli(j+1)/(j+1))
    return re;


############################################################################################
##    Section II: Expanding a given polynomial in terms of a certain basis                ##
############################################################################################

def MonicExpPoly(inPoly, variable, BaseList):     
    '''
    Given 

    inPoly = a polynomials, say P;
    variable = the single variable of P, say y;
    BaseList = the list of polynomials, which form and a basis, denoted by L,

    MonicExpPoly(P, y, L) returns the coefficients of the expansion of P, in terms of L. 

    Note that the polynomials in L, MUST be MONIC and listed in the ASCENDING order. 
    This acturally suggests that L[0] = 1. 



    For example, let 

    P = y^3 + 2 * y - 3, 

    and consider the basis

    L = [1, y - 2, y^2 +1, y^3 + y]. 

    Then, MonicExpPoly(P, y, L) = [-1, 1, 0, 1]
    since (-1) *  + (y - 2) + (y^3 + y) = y^3 + 2 * y + 3. 

    '''

    f = inPoly;                            # polynomial f #
    y = variable;                          # varialbe localy #
    d = f.degree(y)                        # degree of the polynomial, which shall be used for the loop #
    re = []                                # return coefficients #
    while d > 0:
        temp = f.coefficient(y, d)         # coefficient of y^d #
        re = [temp] + re                   # insert the coefficient at the beginning #
        f -= temp * BaseList[d]            # substract the term and reduce the degree #
        f = expand(f)
        d -= 1
    return [f] + re                        # After all the step, f is only the constant coefficient. #



############################################################################################
##    Section III: Orthogonal Polynomials w. r. t. Bernoulli and Euler Numbers            ##
############################################################################################

def OrthBer(n, y):                                                        #The orthogonal polynomials w. r. t. Bernoulli B
    if n == -1:
        re = 0;
    if n == 0:
        re = 1;
    if n > 0:                                                       
        re = (y + 1 / 2) * OrthBer(n - 1, y) + (n - 1)^4 / 4 / (2 * n - 1) / (2 * n - 3) * OrthBer(n - 2, y);   #Recurrence
    return expand(re)

def OrthBerHalf(n, y):                                                        #The orthogonal polynomials w. r. t. Bernoulli (B+1/2), which is "Even"
    if n == -1:
        re = 0;
    if n == 0:
        re = 1;
    if n > 0:
        re = y * OrthBerHalf(n - 1,y) + (n - 1)^4/4/(2 * n - 1) / (2 * n - 3) * OrthBerHalf(n - 2,y);         #Recurrence
    return expand(re)

def OrthEuler(n, y):                                                        #The orthogonal polynomials w. r. t. Euler numbers
    if n == -1:
        re = 0;
    if n == 0:
        re = 1;
    if n > 0:
        re = y * OrthEuler(n - 1,y) + (n - 1)^2 * OrthEuler(n - 2,y);                           #Recurrence
    return expand(re)

def OrthEulerSym(n,y):                                                        #The orthogonal polynomials w. r. t. Euler symbol
    if n == -1:
        re = 0;
    if n == 0:
        re = 1;
    if n > 0:
        re = y * OrthEulerSym(n - 1,y) + (n - 1)^2 / 4 * OrthEulerSym(n - 2,y);                           #Recurrence
    return expand(re)



