# Enzymes

nt_codes = {'A': 'A',
            'C': 'C',
            'T': 'T',
            'G': 'G',
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'B': ['G', 'C', 'T'],
            'H': ['A', 'C', 'T'],
            'D': ['A', 'G', 'T'],
            'V': ['A', 'G', 'C'],
            'N': ['A', 'C', 'T', 'G']}

cas9 = {'guide_start':5, 'guide_length':20, 'pam_start':25,
        'pam':'NGG', 'context_length':30}

cas12a = {'guide_start':9, 'guide_length':23,
          'pam_start':5, 'pam':'TTTV', 'context_length':34}

