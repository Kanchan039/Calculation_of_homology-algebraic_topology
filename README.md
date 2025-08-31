# Algebraic Topology in Python 🌀

This project implements **simplicial homology** using Python and SymPy.  
It can compute **boundary matrices**, **Smith Normal Forms**, and hence **Betti numbers** / **homology groups** of a simplicial complex.

---

## Features
- Build a simplicial complex from a list of simplices.
- Compute boundary operators ∂k.
- Compute Betti numbers βk.
- Display homology groups Hk in a human-readable form (e.g., `H0 = Z`, `H1 = Z`).

---

## Example Usage

### Circle (1-dimensional loop)
```python
simplices = [
    (0,), (1,), (2,),
    (0, 1), (1, 2), (0, 2)
]
circle = SimplicialComplex(simplices)
print(circle.describe_homology(1))
