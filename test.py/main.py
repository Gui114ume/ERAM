import numpy as np

def arnoldi_iteration(A, b, n):
    m = A.shape[0]

    h = np.zeros((n + 1, n))
    Q = np.zeros((m, n + 1))

    q = b / np.linalg.norm(b)
    Q[:, 0] = q

    for k in range(n):
        v = A.dot(q)
        for j in range(k + 1):
            h[j, k] = np.dot(Q[:, j].conj(), v)  # <-- Q needs conjugation!
            v = v - h[j, k] * Q[:, j]

        h[k + 1, k] = np.linalg.norm(v)
        eps = 1e-12
        if h[k + 1, k] > eps:
            q = v / h[k + 1, k]
            Q[:, k + 1] = q
        else:
            return Q, h
    return Q, h



A = [[1,3,3],
    [-2,11,-2],
    [8,-7,6]]
A = np.array(A)

b = [1,1,1]
b = np.array(b)
n = 3

Q, h  = arnoldi_iteration(A,b,n)

print(h)
print(Q)
