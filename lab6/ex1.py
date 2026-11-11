
import random

# ==========================================
# BASIC DNA FUNCTIONS
# ==========================================

random.seed(67)

DNA = "ACGT"


def random_dna(length):
    return "".join(random.choice(DNA) for _ in range(length))


def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return "".join(
        complement[base] for base in reversed(sequence)
    )


# ==========================================
# CREATE A TRANSPOSABLE ELEMENT
# ==========================================

def create_transposon(tir, internal_sequence, tsd):

    right_tir = reverse_complement(tir)

    return (
        tsd
        + tir
        + internal_sequence
        + right_tir
        + tsd
    )


# ==========================================
# GENERATE ARTIFICIAL DNA
# ==========================================

te1 = create_transposon(
    "ACGTTGCAAC",
    random_dna(18),
    "GGTAC"
)

te2 = create_transposon(
    "GATCCGTAAG",
    random_dna(8) + te1 + random_dna(8),
    "TGCAT"
)

te3 = create_transposon(
    "TTACGCGTAA",
    random_dna(8) + te2 + random_dna(8),
    "AACGT"
)

te4 = create_transposon(
    "CCGATGTTCA",
    random_dna(20),
    "CTAGG"
)

dna_sequence = (
    random_dna(20)
    + te3
    + random_dna(25)
    + te4
    + random_dna(20)
)

print("ARTIFICIAL DNA SEQUENCE")
print(dna_sequence)
print("\nDNA length:", len(dna_sequence))


# ==========================================
# DETECT TRANSPOSABLE ELEMENTS
# ==========================================

def detect_transposons(sequence, tir_length=10, tsd_length=5,
                       min_internal=8):

    results = []

    n = len(sequence)

    for start in range(n - tir_length):

        left_tir = sequence[start:start + tir_length]

        expected_right = reverse_complement(left_tir)

        for right_start in range(
            start + tir_length + min_internal,
            n - tir_length + 1
        ):

            right_tir = sequence[
                right_start:right_start + tir_length
            ]

            if right_tir != expected_right:
                continue

            if start < tsd_length:
                continue

            if right_start + tir_length + tsd_length > n:
                continue

            left_tsd = sequence[
                start - tsd_length:start
            ]

            right_tsd = sequence[
                right_start + tir_length:
                right_start + tir_length + tsd_length
            ]

            if left_tsd != right_tsd:
                continue

            element_start = start + 1
            element_end = right_start + tir_length

            results.append({
                "start": element_start,
                "end": element_end,
                "left_TIR": left_tir,
                "right_TIR": right_tir,
                "TSD": left_tsd,
                "length": element_end - element_start + 1
            })

    return results


# ==========================================
# DISPLAY DETECTED ELEMENTS
# ==========================================

detected = detect_transposons(dna_sequence)

print("\nDETECTED TRANSPOSABLE ELEMENTS")
print("-" * 60)

for i, element in enumerate(detected, 1):

    print(f"\nTransposon {i}")
    print("Start position:", element["start"])
    print("End position:", element["end"])
    print("Length:", element["length"])
    print("Left TIR:", element["left_TIR"])
    print("Right TIR:", element["right_TIR"])
    print("TSD:", element["TSD"])


# ==========================================
# DETECT OVERLAPPING ELEMENTS
# ==========================================

print("\nOVERLAPPING TRANSPOSABLE ELEMENTS")
print("-" * 60)

overlap_found = False

for i in range(len(detected)):

    for j in range(i + 1, len(detected)):

        a = detected[i]
        b = detected[j]

        if (a["start"] <= b["end"] and
                b["start"] <= a["end"]):

            overlap_start = max(a["start"], b["start"])
            overlap_end = min(a["end"], b["end"])

            overlap_found = True

            print(
                f"Transposon {i + 1} and "
                f"Transposon {j + 1} overlap."
            )

            print(
                "Overlapping region:",
                overlap_start,
                "-",
                overlap_end
            )

if not overlap_found:
    print("No overlapping elements detected.")
