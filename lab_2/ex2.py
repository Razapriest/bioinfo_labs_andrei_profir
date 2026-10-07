import math

S = input("Enter DNA sequence: ")

Na = float(input("Enter Na+ concentration: "))

length = len(S)
CG = S.count('C') + S.count('G')
CG_percent = (CG / length) * 100

Tm = 81.5 + 16.6 * math.log10(Na) + 41 * (CG_percent / 100) - (600 / length)

print("Tm:", round(Tm, 2))