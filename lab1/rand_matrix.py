import random

matrix = []
n = random.randint(1, 20)
for i in range(n):
    row = [random.randint(1, 10) for _ in range(n+1)]
    row[i] = random.randint(1000, 2000)
    matrix.append(" ".join(map(str, row)))

print("\n".join(matrix))