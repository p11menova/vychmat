import os

global N, EPS


def print_M(M):
    for i in M:
        for j in i:
            print(str(round(j, 4)).center(6), end=" ")
        print()
    print()


def read_from_file():
    with open("matrix.txt") as f:
        n = int(f.readline())
        if n > 20:
            raise ValueError

        eps = float(f.readline())
        M = [i.rstrip().split(" ") for i in f.readlines()]
    return eps, n, M


def validate_matrix(M):
    if len(M) != N:
        raise ValueError

    for i in range(len(M)):
        if (len(M[i])) != N + 1:
            raise ValueError
        M[i] = list(map(float, M[i]))
    return M


#  в матрицe N строк и N+1 столцов
def check_diagonal_dominance(M):
    for i in range(N):
        if not (abs(M[i][i]) >= sum([abs(M[i][j]) for j in range(N) if j != i])):
            return False
    return True


def change_to_get_diagonal_dominance(M):
    X = ['x' + str(i) for i in range(N)]  # чтоб сохранить коэфф при переменных
    for i in range(N):
        for j in range(i, N):
            M[i], M[j] = M[j], M[i]

            # поменяла местами строки и теперь для каждой такой новой матрицы надо поменять местами столбцы

            for k in range(N):
                for h in range(k, N):
                    for _ in range(N):
                        M[_][k], M[_][h] = M[_][h], M[_][k]

                    if check_diagonal_dominance(M):
                        X[k], X[h] = X[h], X[k]
                        return X, M
    return None, None


def count_norm(M):
    res = 0
    for i in M:
        res = max(res, sum([abs(j) for j in i]))

    return res


def get_C_D_matrix(M):
    C = []
    D = []
    for i in range(N):
        a_ii = M[i][i] if M[i][i] != 0 else 1
        c = []
        for j in range(N + 1):
            if j == N:
                D.append(M[i][j] / a_ii)
                continue
            if i == j:
                c.append(0)
            else:
                c.append(-1 * M[i][j] / a_ii)
        C.append(c)
        if all([c1==0 for c1 in c]):
            if D[-1] != 0:
                print("cистема не совместна! решений НЕТ")
                exit(0)
            if D[-1] == 0:
                print("да тут ваще бесконечно количество решений..")
                exit(0)

    return C, D


def read_M():
    start = input("введите f, если хотите прочитать данные из файла matrix.txt и любой другой символ в ином случае:")

    if start.strip() == "f" and os.path.exists("matrix.txt"):
        try:
            eps, n, M = read_from_file()

            return eps, n, M

        except ValueError:
            print("невалидные входные данные...\nзавершение работы программы..")
            exit(1)
    else:
        n = int(input("введите размерность матрицы (≤20): "))
        eps = float(input("введите точность: "))
        print("введите матрицу:")
        M = [input().strip().split() for _ in range(n)]
        return eps, n, M


def get_tochnost(x_k, x_k1):
    # max(abs(x_k_i-x_k1_i))
    a = [abs(x_k[i] - x_k1[i]) for i in range(len(x_k))]
    return max(a)


def put_in_correct_order(x_k, X):
    res = [0] * N
    for i in range(N):
        res[i] = round(x_k[X.index(f'x{i}')], 5)
    return res


def gauss_zeidel(C, D, X):
    print("\nРАСЧЁТЫ ПО МЕТОДУ ГАУССА-ЗЕЙДЕЛЯ:\n")
    k = 0
    # x_k это D
    # x_k+1 который мы сча будем заполнять это

    x_k = D  # имеет вид d1 d2 d3
    x_k1 = [0] * N
    print(f"  номер итерации    |{'x^(k)'.center(32)}|  max(|x^(k)_i-x^(k+1)_i|) ")
    print("-" * 90)
    print(f"      0             |{str(put_in_correct_order(x_k, X)).center(32)}|        -         ")
    while k == 0 or not (get_tochnost(x_k, x_k1) <= EPS):
        x_k, x_k1 = x_k1, [0] * N
        # print(f"x({k})", x_k)
        for i in range(N):
            x = 0
            # всего N шагов и на i-том мы берем i-1 слагаемых из x_k+1
            for j in range(i, N):
                x += C[i][j] * x_k[j]

            for h in range(0, i):
                x += C[i][h] * x_k1[h]

            x += D[i]
            x_k1[i] = x
        # здесь получили x_k1
        k += 1
        print(
            f"      {k}             |{str(put_in_correct_order(x_k1, X)).center(32)}|        {round(get_tochnost(x_k, x_k1), 4)}         ")

    return k, x_k1


def main():
    global N, EPS
    try:
        EPS, N, M = read_M()
        M = validate_matrix(M)
        print("БЫЛО")
        print(''.join([('x' + str(i)).center(6) for i in range(N)]))
        print_M(M)
        X, diagonalized_M = change_to_get_diagonal_dominance(M)

        if not diagonalized_M:
            print("невозможно добиться диагонального преобладания(..\n")
            print("проверяем значение нормы матрицы", count_norm(M))
            if count_norm(M) >= 1:
                print("ну вот( условие сходимости не выполнено(\nсистема не сойдется к ответу данным численным методом "
                      "или будет делать это ОЧЕНЬ долго(")
                exit(1)
            print("ура! условие сходимости выполнено! продолжаем..")
            print()
        if diagonalized_M:
            print("СТАЛО для диагонального преобладания")
            print(''.join([x.center(6) for x in X]))
            print_M(diagonalized_M)

        C, D = get_C_D_matrix(diagonalized_M)
        print("получаем матрицу С")
        print_M(C)
        print("и столбик свободных коэфф D:\n" + '\n'.join(list(map(str, D))))

        k, x_k1 = gauss_zeidel(C, D, X)

        print("ИТОГО:")
        x_k1 = put_in_correct_order(x_k1, X)

        print(
            f"\t-потребовалось {k} итераций \n\t-ответ с учетом погрешности: {x_k1}\n\t-округленный ответ: {list([round(i, 3) for i in x_k1])}")

    except ValueError:
        print("невалидные входные данные...\nзавершение работы программы..")
        exit(1)


main()
