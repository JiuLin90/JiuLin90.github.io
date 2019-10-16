def EvalB(expression,variableList,degreeList):
    re=expression;
    n=len(variableList);
    for i in range(n):
        for j in range(degreeList[i])[::-1]:
	    re=re.subs(variableList[i]^(j+1)==bernoulli(j+1));
    return re;


def EvalE(expression,variableList,degreeList):
    re=expression;
    n=len(variableList);
    for i in range(n):
        for j in range(degreeList[i])[::-1]:
	    re=re.subs(variableList[i]^(j+1)==euler_number(j+1));
    return re;

