#!/usr/bin/env sage
from sage.all import *


def est(n, logq):
    alpha = 0.05
    beta = 0.33
    gamma = 17.88
    delta = 0.65

    return -log(alpha * logq / n, 2) * beta * n / logq + gamma * pow(logq / n, delta) * log(n / logq, 2)


def genprime(init, n, reverse=True):
    if reverse:
        p = init - (init % (2 * n)) - 2 * n + 1
        assert p < init and p % (2 * n) == 1
    else:
        p = init - (init % (2 * n)) + 2 * n + 1
        assert p > init and p % (2 * n) == 1

    while not is_prime(p):
        if reverse:
            p -= 2 * n
        else:
            p += 2 * n

    return p


N = (1 << 13)
t = genprime(1 << 16, N, False)
b = 50
l = 3
w = 3

q = [genprime(1 << b, N)]
while len(q) < l:
    q.append(genprime(q[-1], N))

k = int(ceil(l / w))
P = [genprime(q[-1], N)]
while len(P) < k:
    P.append(genprime(P[-1], N))

print(f"est: {est(N, (l + k) * b).n():.2f}")
print(f"N:   {N}")
print(f"t:   {t}")
print(f"b:   {b}")
print(f"l:   {l}")
print(f"k:   {k}")
print(f"w:   {w}")

print()
print(f'set(TIIFHE_N "{N}")')
print(f'set(TIIFHE_T "\\"{t}\\"")')
print(f'set(TIIFHE_Q "{",".join(map(str, q))}")')
print(f'set(TIIFHE_P "{",".join(map(str, P))}")')
print(f'set(TIIFHE_OMEGA "{w}")')
