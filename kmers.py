
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    list_kmer = []
    kmer = 0
    reverse_kmer = 0
    
    for i in range(k) :
        kmer = kmer << 2
        kmer += encode(text[i])
        reverse_kmer += (3 - encode(text[i])) * (2**(2*i))
    list_kmer.append(min(kmer, reverse_kmer))

    mask = (1<<(k-1)*2)-1
    for n in text[k:] :
        kmer = kmer & mask
        kmer = kmer << 2
        kmer = kmer + encode(n)
                         
        reverse_kmer = reverse_kmer >> 2
        reverse_kmer += (3 - encode(n)) * (2**(2*k))
        list_kmer.append(min(kmer, reverse_kmer))
        
    return list_kmer


def encode(x) :
    dico = {'A':0, 'C':1, 'T':2, 'G':3}
    if x in dico.keys() :
        return dico[x]
    else :
        return 0
