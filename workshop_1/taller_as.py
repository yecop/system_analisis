from io import DEFAULT_BUFFER_SIZE
import math
from collections import Counter
import random

def create_sequence():
    nucleotid_bases = ["A", "C", "G", "T"]
    size_sequence = random.randint(10, 20)
    new_sequence = [nucleotid_bases[random.randint(0, 3)] for i in range(size_sequence)]
    return "".join(new_sequence)

def create_database():
    db_size = 50
    data_base = [create_sequence() for i in range(db_size)]
    return data_base

def shannon_entropy(sequence):
    length = len(sequence)
    freqs = Counter(sequence)
    entropy = 0.0
    for key in freqs:
        prob = freqs[key] / length
        entropy -= prob * math.log2(prob)
    return entropy

def filter_sequences(sequences, threshold=1.5):
    filtered_sequences = []
    for seq in sequences:
        entropy = shannon_entropy(seq)
        if entropy > threshold:
            filtered_sequences.append(seq)
    return filtered_sequences

def get_combinations(n, sequences, bases):
    if n == 1:
        return [sequence + base for sequence in sequences for base in bases]
    else:
        sequence_ = [sequence + base for sequence in sequences for base in bases]
        return get_combinations(n - 1, sequence_, bases)

def count_motif(motif, sequences_db):
    count = 0
    for sequence in sequences_db:
        count += sequence.count(motif)
    return count

def get_motif(motif_size, sequences_db):
    nucleotid_bases = ['A', 'C', 'G', 'T']
    combinations = get_combinations(motif_size, [""], nucleotid_bases)
    max_counter = 0
    motif_winner = ""
    for motif_candidate in combinations:
        temp_counter = count_motif(motif_candidate, sequences_db)
        if temp_counter > max_counter:
            max_counter = temp_counter
            motif_winner = motif_candidate
    return motif_winner, max_counter

def main():
    sequences_db = create_database()
    print("Base de datos original:", sequences_db)

    # Filtrar secuencias usando Entropía de Shannon
    filtered_sequences = filter_sequences(sequences_db)
    print("Base de datos filtrada:", filtered_sequences)

    # Obtener motivos de tamaño 6 y 8
    motif_size_6 = get_motif(6, filtered_sequences)
    motif_size_8 = get_motif(8, filtered_sequences)

    print("Motif de tamaño 6:", motif_size_6)
    print("Motif de tamaño 8:", motif_size_8)

if __name__ == "__main__":
    main()

