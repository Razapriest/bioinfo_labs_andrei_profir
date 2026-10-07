S = input("Enter DNA sequence: ")

G = S.count('G')
C = S.count('C')
T = S.count('T')
A = S.count('A')

Tm = 4 * (G + C) + 2 * (A + T)

print("Melting temperature =", Tm)