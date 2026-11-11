import matplotlib.pyplot as plt
from collections import Counter
import os

def read_dna(filename):
    with open(filename, "r") as file:
        sequence = file.read().strip().upper()

    return sequence

def reverse_complement(sequence):
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]

def read_enzymes():
    enzymes = []

    n = int(input("Number of restriction enzymes: "))

    for i in range(n):
        print(f"\nEnzyme {i + 1}")

        name = input("Name: ")
        motif = input("Recognition sequence: ").upper()
        cut = int(input("Cut position: "))

        enzymes.append({
            "name": name,
            "motif": motif,
            "cut": cut
        })

    return enzymes


def find_restriction_sites(dna, motif, cut_position):
    cleavage_positions = set()

    for pattern, offset in [
        (motif, cut_position),
        (reverse_complement(motif), len(motif) - cut_position)
    ]:
        start = 0

        while True:
            pos = dna.find(pattern, start)

            if pos == -1:
                break

            cleavage_positions.add(pos + offset)
            start = pos + 1

    return sorted(
        pos for pos in cleavage_positions
        if 0 < pos < len(dna)
    )

def calculate_fragments(dna_length, cleavage_positions):

    boundaries = [0] + sorted(set(cleavage_positions)) + [dna_length]

    fragments = []

    for i in range(len(boundaries) - 1):
        fragment_length = boundaries[i + 1] - boundaries[i]

        if fragment_length > 0:
            fragments.append(fragment_length)

    return fragments

def display_results(enzyme_results):

    print("\n" + "=" * 55)
    print("RESTRICTION DIGESTION RESULTS")
    print("=" * 55)

    for result in enzyme_results:

        print(f"\nEnzyme: {result['name']}")
        print(f"Recognition sequence: {result['motif']}")
        print(f"Cut position: {result['cut']}")

        positions = result["positions"]
        fragments = result["fragments"]

        print(f"Number of cleavages: {len(positions)}")

        print("Cleavage positions:")

        if positions:
            print(", ".join(map(str, positions)))
        else:
            print("No cleavage sites found.")

        print("Fragment lengths:")

        for i, fragment in enumerate(fragments, start=1):
            print(f"Fragment {i}: {fragment}")

        print(f"Total fragments: {len(fragments)}")


def simulate_gel(dna_length, enzyme_results):

    lanes = [("Undigested", [dna_length])]

    for result in enzyme_results:
        lanes.append((result["name"], result["fragments"]))

    all_positions = []

    for result in enzyme_results:
        all_positions.extend(result["positions"])

    combined_fragments = calculate_fragments(
        dna_length, sorted(set(all_positions))
    )

    lanes.append(("Combined", combined_fragments))

    fig, ax = plt.subplots(figsize=(10, 7))

    for i, (name, fragments) in enumerate(lanes):

        for size in set(fragments):
            ax.hlines(size, i - 0.3, i + 0.3, color="black", linewidth=2)

    ax.set_yscale("log")
    ax.set_ylim(dna_length * 1.5, 1)

    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels([name for name, _ in lanes], rotation=30)

    ax.set_ylabel("Fragment size")
    ax.set_title("Electrophoresis Gel Simulation")

    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.show()

def main():

    print("=" * 55)
    print("DNA RESTRICTION DIGESTION SIMULATOR")
    print("=" * 55)

    filename = "dna.txt"

    try:
        dna = read_dna(filename)

        print(f"\nDNA sequence successfully loaded.")
        print(f"DNA length: {len(dna)}")

        enzymes = read_enzymes()

        enzyme_results = []

        for enzyme in enzymes:

            positions = find_restriction_sites(
                dna,
                enzyme["motif"],
                enzyme["cut"]
            )

            fragments = calculate_fragments(
                len(dna),
                positions
            )

            enzyme_results.append({
                "name": enzyme["name"],
                "motif": enzyme["motif"],
                "cut": enzyme["cut"],
                "positions": positions,
                "fragments": fragments
            })

        display_results(enzyme_results)

        print("\nGenerating electrophoresis gel...")

        simulate_gel(len(dna), enzyme_results)

    except (FileNotFoundError, ValueError) as error:
        print(f"\nError: {error}")


if __name__ == "__main__":
    main()
