import matplotlib.pyplot as plt
import csv
import numpy as np

def wczytaj_dane(plik):
    odleglosci = []
    wysokosci = []
    with open(plik, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            odleglosci.append(float(row[0]))
            wysokosci.append(float(row[1]))
    return odleglosci, wysokosci

def interpolacja_lagrange(odleglosci, wysokosci, punkty):
    result = []
    for p in range(len(odleglosci)):
        base = []
        for i in range(len(punkty)):
            basis = 1
            for j in range(len(punkty)):
                if i != j:
                    basis *= (odleglosci[p] - odleglosci[punkty[j]]) / (odleglosci[punkty[i]] - odleglosci[punkty[j]])
            base.append(basis)

        y = 0
        for k in range(len(punkty)):
            y += wysokosci[punkty[k]] * base[k]
        result.append(y)
    return result

def interpolacja_sklejana(odleglosci, wysokosci, punkty):
    n = len(punkty) - 1
    size = 4 * n

    A = [[0] * size for _ in range(size)]
    B = [0] * size

    row = 0

    for i in range(n):
        h = odleglosci[punkty[i + 1]] - odleglosci[punkty[i]]
        A[row][4 * i] = 1
        B[row] = wysokosci[punkty[i]]
        row += 1

        A[row][4 * i] = 1
        A[row][4 * i + 1] = h
        A[row][4 * i + 2] = h ** 2
        A[row][4 * i + 3] = h ** 3
        B[row] = wysokosci[punkty[i + 1]]
        row += 1

    for i in range(1, n):
        h_prev = odleglosci[punkty[i]] - odleglosci[punkty[i - 1]]

        A[row][4 * (i - 1) + 1] = 1
        A[row][4 * (i - 1) + 2] = 2 * h_prev
        A[row][4 * (i - 1) + 3] = 3 * h_prev ** 2
        A[row][4 * i + 1] = -1
        B[row] = 0
        row += 1

        A[row][4 * (i - 1) + 2] = 2
        A[row][4 * (i - 1) + 3] = 6 * h_prev
        A[row][4 * i + 2] = -2
        B[row] = 0
        row += 1

    A[row][2] = 1
    B[row] = 0
    row += 1

    A[row][size - 2] = 2
    A[row][size - 1] = 6 * (odleglosci[punkty[-1]] - odleglosci[punkty[-2]])
    B[row] = 0

    w = np.linalg.solve(A, B)

    result = []

    p = 0

    for i in range(len(odleglosci)):
        if(i > punkty[p+1]):
            p += 1
        tmp = 0
        tmp += w[p*4]
        tmp += w[p*4 + 1] * (odleglosci[i] - odleglosci[punkty[p]])
        tmp += w[p*4 + 2] * (odleglosci[i] - odleglosci[punkty[p]])**2
        tmp += w[p*4 + 3] * (odleglosci[i] - odleglosci[punkty[p]])**3
        result.append(tmp)
    return result

def wykres_l1(odleglosci, wysokosci, interp_l, punkty):
    plt.figure(figsize=(12, 6))

    plt.plot(odleglosci, wysokosci, 'o', label='Dane oryginalne')

    plt.plot(odleglosci, interp_l, label='Interpolacja Lagrange’a')

    interpolowane_odleglosci = [odleglosci[i] for i in punkty]
    interpolowane_wysokosci = [wysokosci[i] for i in punkty]
    plt.plot(interpolowane_odleglosci, interpolowane_wysokosci, 'ro', label='Punkty interpolacji równomiernie rozłożone')

    plt.xlabel('Odległość')
    plt.ylabel('Wysokość')
    plt.legend()
    plt.title('Interpolacja profilu wysokościowego')
    plt.show()

def wykres_l2(odleglosci, wysokosci, interp_l2, punkty2):
    plt.figure(figsize=(12, 6))

    plt.plot(odleglosci, wysokosci, 'o', label='Dane oryginalne')

    plt.plot(odleglosci, interp_l2, label='Interpolacja Lagrange’a')

    interpolowane_odleglosci = [odleglosci[i] for i in punkty2]
    interpolowane_wysokosci = [wysokosci[i] for i in punkty2]
    plt.plot(interpolowane_odleglosci, interpolowane_wysokosci, 'ro',
             label='Punkty interpolacji nierównomiernie rozłożone')

    plt.xlabel('Odległość')
    plt.ylabel('Wysokość')
    plt.legend()
    plt.title('Interpolacja profilu wysokościowego')
    plt.show()

def wykres_s1(odleglosci, wysokosci, interp_s, punkty):
    plt.figure(figsize=(12, 6))

    plt.plot(odleglosci, wysokosci, 'o', label='Dane oryginalne')

    plt.plot(odleglosci, interp_s, label='Interpolacja funkcjami sklejanymi')

    interpolowane_odleglosci = [odleglosci[i] for i in punkty]
    interpolowane_wysokosci = [wysokosci[i] for i in punkty]
    plt.plot(interpolowane_odleglosci, interpolowane_wysokosci, 'ro',
             label='Punkty interpolacji równomiernie rozłożone')

    plt.xlabel('Odległość')
    plt.ylabel('Wysokość')
    plt.legend()
    plt.title('Interpolacja profilu wysokościowego')
    plt.show()


def wykres_s2(odleglosci, wysokosci, interp_s2, punkty2):
    plt.figure(figsize=(12, 6))

    plt.plot(odleglosci, wysokosci, 'o', label='Dane oryginalne')

    plt.plot(odleglosci, interp_s2, label='Interpolacja funkcjami sklejanymi')


    interpolowane_odleglosci = [odleglosci[i] for i in punkty2]
    interpolowane_wysokosci = [wysokosci[i] for i in punkty2]
    plt.plot(interpolowane_odleglosci, interpolowane_wysokosci, 'ro', label='Punkty interpolacji nierównomiernie rozłożone')

    plt.xlabel('Odległość')
    plt.ylabel('Wysokość')
    plt.legend()
    plt.title('Interpolacja profilu wysokościowego')
    plt.show()

plik = '2018_paths\\WielkiKanionKolorado.csv'
odleglosci, wysokosci = wczytaj_dane(plik)

liczba_punktow = 50


punkty= np.linspace(0, len(odleglosci)-1, liczba_punktow, dtype=int)

i = np.arange(1, liczba_punktow - 1)
czebyszew = np.cos((2 * i - 1) / (2 * (liczba_punktow - 2)) * np.pi)
dodatkowe_punkty = ((czebyszew + 1) * (len(odleglosci) - 1) / 2).astype(int)

punkty2 = np.unique(np.sort(np.concatenate(([0], dodatkowe_punkty, [len(odleglosci) - 1]))))

interp_l = interpolacja_lagrange(odleglosci, wysokosci, punkty)
interp_l2 = interpolacja_lagrange(odleglosci, wysokosci, punkty2)
interp_s = interpolacja_sklejana(odleglosci, wysokosci, punkty)
interp_s2 = interpolacja_sklejana(odleglosci, wysokosci, punkty2)

wykres_l1(odleglosci, wysokosci, interp_l, punkty)
wykres_l2(odleglosci, wysokosci, interp_l2, punkty2)
wykres_s1(odleglosci, wysokosci, interp_s, punkty)
wykres_s2(odleglosci, wysokosci, interp_s2, punkty2)
