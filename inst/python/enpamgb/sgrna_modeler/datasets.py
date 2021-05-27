"""Training and testing datasets for modeling guide activity.

Includes class for generating new dataset objects, so new data can be easily modeled and tested

    :Example:

    >>> import pandas as pd
    >>> import sgrna_modeler.enzyme as en
    >>> new_data = pd.read_csv('new dataset')
    >>> new_dataset = Activity_Data(new_data, en.cas9, '30mer', 'activity', 'new data')
"""
import pandas as pd
from sgrna_modeler import enzymes as en
import os


def curr_path():
    return os.path.dirname(__file__)


class ActivityData(object):
    """Store information about activity data


    :param data: data to model
    :type data: pandas dataframe
    :param enzyme: cas9 or cas12a
    :type enzyme: dict
    :param kmer_column: sequences to model
    :type kmer_column: str
    :param name: name of the dataset
    :type name: str
    :param group_column: column to include in prediction output
    :type group_column:str
    """

    def __init__(self, data, enzyme, kmer_column, activity_column, name, group_column = ''):
        """Inits Activity data"""
        self.data = data
        self.enzyme = enzyme
        self.kmer_column = kmer_column
        self.activity_column = activity_column
        self.name = 'D_' + name
        self.group_column = group_column

    def get_xy(self):
        """Gets modeling matrix (x) and output matrix (y)

        :return two series, x and y
        :rtype pandas series
        """
        x = self.data[self.kmer_column]
        y = self.data[self.activity_column]
        return x, y


# SpCas9 Datasets

def load_doench_2016():
    """Data from:

    Doench, John G., et al. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." \
    Nature biotechnology 34.2 (2016): 184.

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> doench = da.load_doench_2016()
    >>> doench.data
          Unnamed: 0                           30mer  ...       drug  predictions
    0              0  CAGAAAAAAAAACACTGCAACAAGAGGGTA  ...     nodrug     0.544412
    1              1  TTTTAAAAAACCTACCGTAAACTCGGGTCA  ...    PLX_2uM     0.617512
    2              2  TCAGAAAAAGCAGCGTCAGTGGATTGGCCC  ...     nodrug     0.476232
    3              3  AATAAAAAATAGGATTCCCAGCTTTGGAAG  ...    PLX_2uM     0.459882
    4              4  GATGAAAAATATGTAAACAGCATTTGGGAC  ...    PLX_2uM     0.290841
              ...                             ...  ...        ...          ...
    5305        5305  GCACTTTGGTGTGGCTGACTGAGTGGGCCA  ...    PLX_2uM     0.586758
    5306        5306  TTCTTTTGTAAGAACCCGCTGTGTTGGTTT  ...    PLX_2uM     0.492066
    5307        5307  GCCCTTTGTCATCGTAGGAAGATATGGCTG  ...  AZD_200nM     0.479728
    5308        5308  CAAATTTGTTCTTTAAATGGCTACAGGAGG  ...  AZD_200nM     0.478952
    5309        5309  CAAATTTGTTCTTTAAATGGCTACAGGAGG  ...    PLX_2uM     0.478952
    [5310 rows x 9 columns]

    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Doench_2016.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='30mer',
                              activity_column='score_drug_gene_rank',
                              name='Doench_2016', group_column='Target gene')
    return data_class


def load_meyers_2017_train():
    """Essential genes from GeckoV2 achilles screens:

    Meyers, Robin M., et al. "Computational correction of copy number effect improves specificity of CRISPR–Cas9 \
    essentiality screens in cancer cells." Nature genetics 49.12 (2017): 1779-1784.

    Mean activity is averaged accross screens after Z-scoring by non-essentials

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> meyers_2017_train = da.load_meyers_2017_train()
    >>> meyers_2017_train.data
             Species   Build  Chromosome Number  ... Percent Protein Notes mean_activity
    0      human  GRCh38                1.0  ...           24.87   NaN     -0.230160
    1      human  GRCh38                1.0  ...           23.26   NaN      3.045755
    2      human  GRCh38                1.0  ...           18.60   NaN      1.307097
    3      human  GRCh38                1.0  ...           18.07   NaN     -1.307698
    4      human  GRCh38                1.0  ...           13.95   NaN      1.278670
          ...     ...                ...  ...             ...   ...           ...
    7897   human  GRCh38                NaN  ...           11.78   NaN      1.959897
    7898   human  GRCh38                NaN  ...           13.14   NaN     -0.429659
    7899   human  GRCh38                NaN  ...           16.01   NaN      1.187820
    7900   human  GRCh38                NaN  ...           19.03   NaN      1.573194
    7901   human  GRCh38                NaN  ...           39.62   NaN      2.044455
    [7902 rows x 20 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Meyers_2017_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='sgRNA context sequence',
                              activity_column='mean_activity',
                              name='Meyers_2017_Train', group_column='Gene Symbol')
    return data_class


def load_meyers_2017_test():
    """Essential genes from GeckoV2 achilles screens:

    Meyers, Robin M., et al. "Computational correction of copy number effect improves specificity of CRISPR–Cas9 \
    essentiality screens in cancer cells." Nature genetics 49.12 (2017): 1779-1784.

    Mean activity is averaged accross screens after Z-scoring by non-essentials

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> meyers_2017_test = da.load_meyers_2017_test()
    >>> meyers_2017_test.data
                Species   Build  Chromosome Number  ... Percent Protein Notes mean_activity
    0     human  GRCh38                1.0  ...           22.12   NaN      3.325952
    1     human  GRCh38                1.0  ...            2.30   NaN      2.645421
    2     human  GRCh38                1.0  ...           56.87   NaN      2.040191
    3     human  GRCh38                1.0  ...           40.38   NaN      3.356250
    4     human  GRCh38                1.0  ...           40.11   NaN      1.602670
    ..      ...     ...                ...  ...             ...   ...           ...
    667   human  GRCh38                NaN  ...            8.56   NaN      1.240547
    668   human  GRCh38                NaN  ...            8.95   NaN      1.078080
    669   human  GRCh38                NaN  ...           30.93   NaN     -0.364154
    670   human  GRCh38                NaN  ...           34.24   NaN      2.605412
    671   human  GRCh38                NaN  ...           41.25   NaN      2.620977
    [672 rows x 20 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Meyers_2017_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9, kmer_column='sgRNA context sequence',
                              activity_column='mean_activity',
                              name='Meyers_2017_Test', group_column='Gene Symbol')
    return data_class


def load_kim_2019_train():
    """
    Indel frequencies from:

    Kim, Hui Kwon, et al. "SpCas9 activity prediction by DeepSpCas9, a deep learning–based model with high \
    generalization performance." Science advances 5.11 (2019): eaax9249.

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> kim_2019_train = da.load_kim_2019_train()
    >>> kim_2019_train.data
                           Barcode  ... Background subtracted indel (%)
    0      TTTGACACACACGCACTAG  ...                       24.287805
    1      TTTGACACACACTCGTATG  ...                       69.500438
    2      TTTGACACACACTCTCGTC  ...                       25.994760
    3      TTTGACACACACTCTGCTG  ...                       57.964590
    4      TTTGACACACACTGCATAT  ...                       39.355020
                        ...  ...                             ...
    12827  TTTGTGTGTCTCGTATCAC  ...                       40.853256
    12828  TTTGTGTGTCTCTACACGC  ...                       11.480880
    12829  TTTGTGTGTCTCTCACGTA  ...                       63.861469
    12830  TTTGTGTGTCTCTCTAGTC  ...                       51.650932
    12831  TTTGTGTGTCTCTCTCAGA  ...                       40.019124
    [12832 rows x 9 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2019_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9,
                              kmer_column='Target context sequence (4+20+3+3)',
                              activity_column='Background subtracted indel (%)',
                              name='Kim_2019_Train')
    return data_class


def load_kim_2019_test():
    """
    Indel frequencies from:

    Kim, Hui Kwon, et al. "SpCas9 activity prediction by DeepSpCas9, a deep learning–based model with high \
    generalization performance." Science advances 5.11 (2019): eaax9249.

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> kim_2019_test = da.load_kim_2019_test()
    >>> kim_2019_test.data
        Target context sequence (4+20+3+3)  ...  Background subtracted indel frequencies (average, %)
    0       AAAACTGTGAGTGTGGGACCTGCTGGGGGC  ...                                          44.125755
    1       AAACACAACCAATCCGAGGCCTTCTGGGTC  ...                                          12.163189
    2       AAACTGTGAGTGTGGGACCTGCTGGGGGCT  ...                                          68.901263
    3       AAACTTGAGAGCTTTCATAAAGCTTGGCAA  ...                                          13.135690
    4       AAAGAAGCGGACTTTAAAGTTCGAGGGAGA  ...                                          48.355156
    ..                                 ...  ...                                                ...
    537     TTTGCAGCGCGTTGACTTATTCATGGGTCA  ...                                          36.249050
    538     TTTGCTAGGAATATTGAAGGGGGCAGGGGA  ...                                          38.622947
    539     TTTGTGGTGGTTGCTATGGTAATCCGGCAC  ...                                          12.246218
    540     TTTTTACAATTCTGTGAGTTAGAGTGGGCA  ...                                           0.385915
    541     TTTTTGAGGTGCACTAATAGAGGGTGGAGT  ...                                          41.100730
    [542 rows x 5 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2019_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas9,
                              kmer_column='Target context sequence (4+20+3+3)',
                              activity_column='Background subtracted indel frequencies\r(average, %)',
                              name='Kim_2019_Test')
    return data_class


# AsCas12a datasets
def load_kim_2018_train():
    """
    Indel frequencies from:

    Kim, Hui Kwon, et al. "Deep learning improves prediction of CRISPR–Cpf1 guide RNA activity." \
    Nature biotechnology 36.3 (2018): 239.

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> kim_2018_train = da.load_kim_2018_train()
    >>> kim_2018_train.data
              50 bp synthetic target and target context sequence 10 bp + PAM + 23 bp protospacer + 17 bp)  ... Indel frequency
    0      TGCGCGAGCGTTTAAAAAACATCGAACGCATCTGCTGCCTAGCTTG...                                             ...       14.711302
    1      CTAAAGAAACTTTAAAAATCTTTTCTGCCAGATCTCCAGAAGCTTG...                                             ...        0.238095
    2      TTGCCATTGTTTTAAAACAGGTTCTGTACTTGATCTCTCCAGCTTG...                                             ...       88.079746
    3      TTGCACATATTTTAAAACTGAGTTCAAAGACCACTCTTCCAGCTTG...                                             ...       75.392670
    4      TAGACTAATGTTTAAAAGCAAGTGCAAGTCTTTGGAATCTAGCTTG...                                             ...       63.320080
                                                      ...                                             ...             ...
    14995  TCCATCTTCATTTTTTTTGTAGAGTAGGGCTTTATTTCCAAGCTTG...                                             ...       -0.467290
    14996  CCTTCTCTCCTTTTTTTTTCAAGATCTGATTCTTCTTGCAAGCTTG...                                             ...        0.000000
    14997  CCAGGACTTGTTTTTTTTTCAATCTGTTCATCTTGGACCAAGCTTG...                                             ...        0.239006
    14998  ACCATCATAATTTTTTTTTGCAACATAGCCATTTCTTTTTAGCTTG...                                             ...       -0.272826
    14999  GAGCGCTTCTTTTTTTTTTTCGGGGTCTCGTTGCTGGGCGAGCTTG...                                             ...       -2.766164
    [15000 rows x 10 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2018_Train.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas12a,
                              kmer_column='Context Sequence',
                              activity_column='Indel frequency',
                              name='Kim_2018_Train')
    return data_class


def load_kim_2018_test():
    """
    Indel frequencies from:

    Kim, Hui Kwon, et al. "Deep learning improves prediction of CRISPR–Cpf1 guide RNA activity." \
    Nature biotechnology 36.3 (2018): 239.

    :Example:

    >>> import sgrna_modeler.datasets as da
    >>> kim_2018_test = da.load_kim_2018_test()
    >>> kim_2018_test.data
             50 bp synthetic target and target context sequence  ... Indel frequency
    0     GCAATTTGGTTTTAAAACAGAATATACAGTCTAAAAAACCAGCTTG...  ...       71.580711
    1     CTGATGGCCATTTAAACAACTCTTTGAGCTCTCCAGTTCAAGCTTG...  ...       19.672949
    2     TTTAGATGATTTTAAACCAGCATCTATAGACACTTCCTGTAGCTTG...  ...       75.641026
    3     ACATTTGGACTTTAAACCCAAACTACTTGTCCAACGGTACAGCTTG...  ...       46.920217
    4     CTCTACCAGGTTTAAACGCTTCCACACTTGTGTCAGTAATAGCTTG...  ...       54.981550
                                                     ...  ...             ...
    2958  AGTTTGGAATTTTTTTTACACTGATCCTCAGCACATCTCAAGCTTG...  ...       -0.378500
    2959  CAGGCTTTCTTTTTTTTCCTTTCCTAGTTGGTTCATTCCCAGCTTG...  ...        0.189438
    2960  AACAGTGGCTTTTTTTTGCTGCTAGCACATATGTATGGGTAGCTTG...  ...       -2.857143
    2961  CAGCCTCATGTTTTTTTGGGAACCAATCGATAATCACATTAGCTTG...  ...       11.275673
    2962  TTGGATTGTGTTTTTTTTTAGCACCTTATTTTCCTTGAAGAGCTTG...  ...       -1.675978
    [2963 rows x 10 columns]
    """
    data = pd.read_csv(os.path.join(curr_path(), 'data/datasets/Kim_2018_Test.csv.zip'))
    data_class = ActivityData(data=data, enzyme=en.cas12a,
                              kmer_column='Context Sequence',
                              activity_column='Indel frequency',
                              name='Kim_2018_Test')
    return data_class

# enAsCas12a datasets
