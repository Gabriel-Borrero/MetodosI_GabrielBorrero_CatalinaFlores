import sympy as sp

# Definir las matrices γμ
gamma0 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
gamma1 = sp.Matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]])
gamma2 = sp.Matrix([[0, 0, 0, -sp.I], [0, 0, sp.I, 0], [0, sp.I, 0, 0], [-sp.I, 0, 0, 0]])
gamma3 = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]])

# Calcular {γμ, γν}
gamma_mu = [gamma0, gamma1, gamma2, gamma3]
anti_commutation = []
for mu in range(4):
    for nu in range(4):
        result = gamma_mu[mu] * gamma_mu[nu] + gamma_mu[nu] * gamma_mu[mu]
        anti_commutation.append(result)

# Definir la métrica ημν
eta = sp.diag(1, -1, -1, -1)

# Verificar la relación de anticonmutación
for mu in range(4):
    for nu in range(4):
        expected_result = 2 * eta[mu, nu] * sp.eye(4)
        print(expected_result)
        if anti_commutation[mu * 4 + nu] == expected_result:
            print(f"{{γ{mu}, γ{nu}}} = {expected_result} (Verificado)")
        else:
            print(f"{{γ{mu}, γ{nu}}} ≠ {expected_result} (No verificado)")
