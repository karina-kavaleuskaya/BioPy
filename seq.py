from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations
import random

# Set the length of the random DNA sequence
sequence_length = 10000

# Create a random DNA sequence
dna_sequence = ''.join(random.choices('ATGC', k=sequence_length))

# Create a Seq object from the DNA sequence
seq = Seq(dna_sequence)
print(seq)

# Get the six frame translations for the DNA sequence
#translations = six_frame_translations(str(seq))

# Create a dictionary to store the translations
#translations_dict = {}

# Populate the dictionary with frame numbers and translations
#for i, translation in enumerate(translations):
#    frame = i+1
#    translations_dict[frame] = str(translation)

# Print the protein sequences for all frame translations
#for frame, translation in translations_dict.items():
#    print(f"Frame {frame}: {translation}")