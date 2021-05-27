import pandas as pd
import numpy as np
import random

def generate_variants(seq):
  # generate a list of all possible variants of a sequence
  variants = []
  variant_nt = []
  variant_pos = []
  nts = ['A', 'C', 'T', 'G']
  for i, seq_nt in enumerate(seq):
    for N in nts:
      if seq_nt != N:
        new_seq = seq[:i] + N + seq[i+1:]
        variants.append(new_seq)
        variant_nt.append(N)
        variant_pos.append(i)
  variant_df = pd.DataFrame({'variant': variants, 'nt': variant_nt, 'pos': variant_pos})
  return variant_df

def mutagenize_seq(guide, model):
  variant_df = generate_variants(guide)
  variant_predictions = model.predict_seqs(variant_df.variant)
  variant_df['prediction'] = variant_predictions
  original_prediction = model.predict_seqs([guide])
  variant_df['delta'] = variant_df.prediction - original_prediction
  summarized_df = variant_df.groupby('pos').agg({'delta': 'mean'}).reset_index()
  summarized_df['nt'] = list(guide)
  summarized_df['importance'] = np.absolute(summarized_df.delta)
  summarized_df['pos'] = summarized_df.pos + 1
  summarized_df['context'] = guide
  return summarized_df


def mutagenize_model(model, resamples):
  nts = ['A', 'C', 'T', 'G']
  seq_len = model.enzyme['context_length']
  mutation_list = []
  sequence_list = []
  sequence = ''.join(random.choice(nts) for i in range(seq_len))
  mutation = {'sequence': sequence, 'reference': '', 'nt': '','position': float('nan')}
  mutation_list.append(mutation)
  sequence_list.append(sequence)
  mutation_list = []
  for i in range(resamples):
    while sequence in sequence_list:
      last_sequence = sequence
      position = random.randint(0, seq_len - 1)
      nt = random.choice(nts)
      sequence = sequence[:position] + nt + sequence[(position + 1):]
      mutation = {'sequence': sequence, 'reference': last_sequence, 'nt': nt, 'position': position + 1}
    mutation_list.append(mutation)
    sequence_list.append(sequence)
  sequence_df = pd.DataFrame(mutation_list)
  sequence_predictions = model.predict_seqs(sequence_df.sequence)
  prediction_df = pd.DataFrame({'sequence': sequence_df.sequence, 'prediction': sequence_predictions})
  delta_df = (sequence_df.merge(prediction_df, on = 'sequence')
              .merge(prediction_df, left_on = 'reference', right_on = 'sequence', suffixes = ('_seq', '_ref')))
  delta_df['delta'] = delta_df.prediction_seq - delta_df.prediction_ref
  return delta_df

  

