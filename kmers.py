import heapq 


def kmer2str(val, k):
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for i in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def stream_kmers(text, k):
    # Precompute the (k-1)-mer (and reverse)
    kmer = 0
    rkmer = 0
    for letter in text[:k-1]:
        # A = 00, C = 01, T = 10, G = 11
        # Forward kmer
        kmer <<= 2
        letter_value = (ord(letter) >> 1) & 0b11
        kmer += letter_value
        # Reverse kmer
        rkmer >>= 2
        rev_letter_value = (letter_value + 2) & 0b11
        rkmer += rev_letter_value << (2 * (k - 1))

    # Stream kmers
    mask = (1 << (2 * k)) - 1
    for letter in text[k-1:]:
        # Forward kmer
        kmer <<= 2
        letter_value = (ord(letter) >> 1) & 0b11
        kmer += letter_value
        kmer &= mask
        # Reverse kmer
        rkmer >>= 2
        rev_letter_value = (letter_value + 2) & 0b11
        rkmer += rev_letter_value << (2 * (k - 1))

        yield kmer, rkmer



        
def sketch(s, stream, seq, k) :
    L = []
    for kmer, rkmer in stream(seq, k) :
        Kmer = min(kmer, rkmer)
        Kmer = xorshift(Kmer)
        
        if len(L) < s :
            L.append(Kmer)
        else :
            Max = max(L)
            if Kmer > Max :
                index_max = L.index(Max)
                L[index_max] = Kmer
                
    return sorted(L)




def sketch_queue(s, stream, seq, k) :
    L = []
    
    for kmer, rkmer in stream(seq, k) :
        Kmer = min(xorshift(kmer), xorshift(rkmer))
        if len(L) < s :
            L.append(Kmer)
            if len(L) == s :
                heapify(L)
        else :
            Max = maxheap(L)
            if Kmer > Max :
                heappop(L)
                heappush(L, Kmer)
                
    return sorted(L)
            

    
def xorshift(x) :
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return x


def hash(list_kmers) :
    return [xorshift(i) for i in list_kmers]


print(sketch(100, stream_kmers, 'CGGAATTA', 3))
            
