import sympy as sp
from sympy.matrices.normalforms import smith_normal_form as _snf
from sympy.polys.domains import ZZ

class SimplicialComplex:
    def __init__(self, simplices):
        self.simplices = simplices
        self.max_dim = max(len(s) for s in simplices) - 1
        self.simplices_by_dim = {k: [] for k in range(self.max_dim + 1)}
        for s in simplices:
            self.simplices_by_dim[len(s) - 1].append(tuple(sorted(s)))

        for k in self.simplices_by_dim:
            self.simplices_by_dim[k] = list(sorted(set(self.simplices_by_dim[k])))

    def boundary_matrix(self, k):
        """Return boundary matrix ∂_k: C_k -> C_{k-1}"""
        if k <= 0 or k > self.max_dim:
            return sp.zeros(0)

        cols = self.simplices_by_dim[k]
        rows = self.simplices_by_dim[k - 1]

        M = sp.zeros(len(rows), len(cols))

        for j, sigma in enumerate(cols):
            for i, v in enumerate(sigma):
                face = tuple(sorted(sigma[:i] + sigma[i+1:]))
                sign = (-1) ** i
                if face in rows:
                    M[rows.index(face), j] = sign
        return M

    def smith_normal_form(self, A: sp.Matrix):
        """Compute Smith Normal Form over Z (works across SymPy versions)."""
        if A.rows == 0 or A.cols == 0:
            D = sp.zeros(A.rows, A.cols)
            return sp.eye(A.rows), D, sp.eye(A.cols)

        try:
            D, U, V = _snf(A, domain=ZZ, calc_transformation=True)
            return U, D, V
        except TypeError:
            D = _snf(A, domain=ZZ)
            U = sp.eye(A.rows)
            V = sp.eye(A.cols)
            return U, D, V

    def homology(self, k):
        """Compute Betti number of H_k"""
        Bk = self.boundary_matrix(k)
        Bk1 = self.boundary_matrix(k + 1)

        # Rank of Bk (image of ∂_k)
        _, Dk, _ = self.smith_normal_form(Bk)
        rank_Bk = sum(1 for i in range(min(Dk.rows, Dk.cols)) if Dk[i, i] != 0)

        # Rank of Bk+1 (image of ∂_{k+1})
        _, Dk1, _ = self.smith_normal_form(Bk1)
        rank_Bk1 = sum(1 for i in range(min(Dk1.rows, Dk1.cols)) if Dk1[i, i] != 0)

        dimCk = len(self.simplices_by_dim.get(k, []))
        betti_k = dimCk - rank_Bk - rank_Bk1
        return betti_k

    def betti_numbers(self, up_to=None):
        """Return Betti numbers β_0 ... β_up_to"""
        if up_to is None:
            up_to = self.max_dim
        return [self.homology(k) for k in range(up_to + 1)]

    def describe_homology(self, up_to=None):
        """Return human-readable homology groups"""
        bettis = self.betti_numbers(up_to)
        result = []
        for k, b in enumerate(bettis):
            if b == 0:
                group = "0"
            elif b == 1:
                group = "Z"
            else:
                group = "Z^{}".format(b)
            result.append(f"H{k} = {group}")
        return "\n".join(result)


# ---------------- Example usage ----------------
if __name__ == "__main__":
    # Circle (1-dimensional simplicial complex: 3 vertices, 3 edges)
    simplices = [
        (0,), (1,), (2,),
        (0, 1), (1, 2), (0, 2)
    ]
    circle = SimplicialComplex(simplices)
    print("Circle Complex:")
    print(circle.describe_homology(1))
    print()

    # Filled triangle (2-simplex: contractible)
    simplices2 = [
        (0,), (1,), (2,),
        (0, 1), (1, 2), (0, 2),
        (0, 1, 2)
    ]
    triangle = SimplicialComplex(simplices2)
    print("Triangle Complex:")
    print(triangle.describe_homology(2))
