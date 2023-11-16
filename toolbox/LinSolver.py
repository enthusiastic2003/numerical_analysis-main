def printMatrix(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            print(str(A[i][j])+" ",end="")
        print("\n")

def GaussElem(A,B):
    #Convert to augmentd matrix
    AugM=[]
    n=len(B)
    printMatrix(A)
    AugM = [row + [bi] for row, bi in zip(A, B)]

    #First select a pivot and pivot element and swap rows
    for i in range(n):
        pivot_row=max(range(i,n),key=lambda l: abs(AugM[l][i]))
        AugM[pivot_row],AugM[i]=AugM[i],AugM[pivot_row]
        pivot_element=AugM[i][i]

        AugM[i]=[x/pivot_element for x in AugM[i]]
        #Elemination
        for j in range(i+1,n):
            fact=AugM[j][i]
            AugM[j]=[x - fact * y for x, y in zip(AugM[j], AugM[i])]
    w=[None for i in range(n)]
    for i in range(n-1,-1,-1):
        sum=AugM[i][-1]
        for j in range(n-1,i,-1):
            sum=sum-AugM[i][j]*w[j]
        w[i]=sum/AugM[i][i]
    return w

def GaussJacobi(A, B, N, Tol, X0):
    """
    This Method Solves A system of linear equations given in A.w=B format, using the GaussJacobi Formulation.
    Given: 
    A.w=B
    dim(A)=n X n
    N=How many iterations?
    X0=initial approximation
    """
    n=len(A)
    k = 1
    while k <= N:
        X = [0] * len(A)
        for i in range(n):
            Sum = 0.0
            for j in range(n):
                if j != i:
                    Sum += A[i][j] * X0[j]
            X[i] = (B[i][0] - Sum) / A[i][i]

        print(k, 'th iteration', X)
        diff = max(abs(X[i] - X0[i]) for i in range(len(X)))
        if diff < Tol:
            break
        k += 1
        X0 = X[:]

    if k == N + 1:
        print("The method failed after N0 iterations")
    else:
        print("The method is successful")


def Gauss_Seidal(A, B, N, Tol, X0):
    n=len(A)
    k = 1
    while k <= N:
        X = [0] * len(A)  # Initialize X
        for i in range(n):
            Sum1 = sum(A[i][j] * X[j] for j in range(i))
            Sum2 = sum(A[i][j] * X0[j] for j in range(i + 1, n))
            X[i] = (B[i][0] - Sum1 - Sum2) / A[i][i]

        print(k, 'th iteration', X)
        diff = max(abs(X[i] - X0[i]) for i in range(len(X)))
        if diff < Tol:
            break
        k += 1
        X0 = X[:]

    if k == N + 1:
        print("The method failed after N iterations")
    else:
        print("The method is successful")