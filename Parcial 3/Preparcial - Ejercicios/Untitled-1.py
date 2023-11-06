import numpy as np

# Definir las matrices de Dirac y la métrica del espacio
gamma0 = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, -1, 0],
                  [0, 0, 0, -1]])

gamma1 = np.array([[0, 0, 0, 1],
                  [0, 0, 1, 0],
                  [0, -1, 0, 0],
                  [-1, 0, 0, 0]])

gamma2 = np.array([[0, 0, 0, -1j],
                  [0, 0, 1j, 0],
                  [0, 1j, 0, 0],
                  [-1j, 0, 0, 0]])

gamma3 = np.array([[0, 0, 1, 0],
                  [0, 0, 0, -1],
                  [-1, 0, 0, 0],
                  [0, 1, 0, 0]])

# Métrica del espacio
eta = np.diag([1, -1, -1, -1])

matrices_de_dirac = [gamma0, gamma1, gamma2, gamma3]
results = []

for i in range(4):
    for j in range(4):
        matrix1 = matrices_de_dirac[i]
        matrix2 = matrices_de_dirac[j]

        anticommutation_relation = matrix1 @ matrix2 + matrix2 @ matrix1

        is_equal_to_eta = np.allclose(anticommutation_relation, 2 * eta)
        
        results.append((i, j, anticommutation_relation, is_equal_to_eta))

# Mostrar los resultados
for i, j, anticommutation_relation, is_equal_to_eta in results:
    print(f"Para gamma{i} y gamma{j}:")
    print("Matriz resultante de anticonmutación:")
    print(anticommutation_relation)
    print("¿Igual a 2*eta?:", is_equal_to_eta)
    print()
