"""
Lempel-Ziv complexity measures for binary sequences.

Implements two algorithms:
  - LZ76 (exhaustive): Lempel & Ziv (1976) "On the Complexity of Finite Sequences"
  - LZ78 (dictionary): Lempel & Ziv (1978) compression-based dictionary complexity

Translated from MATLAB calc_lz_complexity.m (Quang Thai, 2012) included in CopBET.
"""
import numpy as np


def lz76_complexity(S, normalize=True):
    """
    Compute LZ76 exhaustive Lempel-Ziv complexity of a binary sequence.

    Parameters
    ----------
    S : array-like, bool or int
        Binary sequence (0s and 1s).
    normalize : bool
        If True, normalize by L / log2(L).

    Returns
    -------
    C : float
        LZ76 complexity value.
    """
    S = np.asarray(S, dtype=bool)
    L = len(S)
    if L == 0:
        return 0.0

    # Build string representation for fast substring search
    S_str = ''.join('1' if x else '0' for x in S)

    # Compute eigenfunction gs[n]: gs[n] = max m such that S[m:n] not in vocab(S[0:n-1])
    gs = [0] * (L + 1)
    gs[0] = 0
    gs[1] = 1

    for n in range(2, L+1):
        idx_list = list(range(gs[n-1]+1, n + 1))  # gs[n-1] (0-indexed) to n
        eigenvalue_found = False

        for k in range((len(idx_list)+1) // 2):
            # Check upper end
            m_upper = idx_list[-(k + 1)]
            if S_str[m_upper-1:n] not in S_str[:n-1]:
                gs[n] = m_upper
                eigenvalue_found = True
                break

            # Check lower end
            # if k < len(idx_list) - k - 1:
            m_lower = idx_list[k]
            if S_str[m_lower-1:n] in S_str[:n-1]:
                gs[n] = m_lower - 1
                eigenvalue_found = True
                break
            # elif idx_list[-(k + 1)] == idx_list[k] + 1:
            elif m_upper == m_lower + 1:
                gs[n] = m_lower
                eigenvalue_found = True
                break

        if not eigenvalue_found:
            raise RuntimeError(f"Failed to find eigenvalue for n={n} in LZ76 complexity calculation.")
            # gs[n + 1] = gs[n]

    # Find exhaustive terminal points
    # h_i = [0]
    # h_prev = 0
    # k = 1
    # while True:
    #     found = False
    #     for h in range(h_prev + 2, L + 2):  # +2 because gs is 1-indexed offset
    #         if h <= L and gs[h] > h_prev:
    #             h_prev = h - 1  # convert back to sequence index
    #             h_i.append(h_prev)
    #             found = True
    #             break
    #     if not found:
    #         break
    gs = np.asarray(gs)
    h_i = np.zeros(len(gs), dtype=int)
    h_i_length = 1  # MATLAB: first element already "present"
    h_prev = 0      # corresponds to h_0

    h_i[0] = h_prev

    k = 1    
    while k is not None:
        # MATLAB:
        # k = find(gs((h_prev+1+1):end) > h_prev, 1);

        start = h_prev + 1  # convert MATLAB (h_prev+2)-1 to Python index

        if start >= len(gs):
            break

        sub = gs[start:]

        # find first index where condition holds
        idxs = np.where(sub > h_prev)[0]

        if len(idxs) == 0:
            break

        k = idxs[0] + 1  # MATLAB find returns 1-based index

        h_i_length += 1
        h_prev = h_prev + k
        h_i[h_i_length - 1] = h_prev
    
    # h_i = h_i[:h_i_length]


    if h_i[h_i_length-1] != L:
        # h_i.append(L)
        h_i_length += 1
        h_i[h_i_length - 1] = L

    C_raw = h_i_length - 1  # number of production components

    if normalize:
        if L > 1:
            C = C_raw / (L / np.log2(L))
        else:
            C = float(C_raw)
    else:
        C = float(C_raw)

    return C


def lz78_complexity(S):
    """
    Compute LZ78 dictionary (compression) complexity of a binary sequence.

    This is the 'cpr' algorithm from Schartner's Python code used in the original
    CopBET MATLAB implementation.

    Parameters
    ----------
    S : array-like, bool or int
        Binary sequence.

    Returns
    -------
    C : int
        Number of dictionary entries (raw complexity count).
    """
    S = np.asarray(S, dtype=bool)
    s = ''.join('1' if x else '0' for x in S)

    dictionary = set()
    w = ''
    count = 0

    for c in s:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            dictionary.add(wc)
            count += 1
            w = c

    return count
