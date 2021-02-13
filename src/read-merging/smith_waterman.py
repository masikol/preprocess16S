# -*- coding: utf-8 -*-
#      |===== Smith-Waterman algorithm implementation =====|

from array import array
from functools import partial

class SWAlignLensNotEq(Exception):
    # An exception meant to be raised when length of the first aligned 
    # sequence the second's one are not equal after aligning.
    # Is subclass of Exception.
    pass
# end class SWAlignLensNotEq


class SWInvalidMove(Exception):
    # An exception meant to be raised when Smith-Waterman algorithm can't perform the next move.
    # Is subclass of Exception.
    pass
# end class SWInvalidMove


class AlignResult:
    # Class AlignResult is dedicated to perform result of Smith-Waterman aligning.
    # :field q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
    # :type q_align: str;
    # :field s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
    # :type s_align: str;
    # :field q_start: 1-based number of query sequence, in which alignment starts;
    # :type q_start: int;
    # :field q_end: 1-based number of query sequence, in which alignment ends;
    # :type q_end: int;
    # :field s_start: 1-based number of subject sequence, in which alignment starts;
    # :type s_start: int;
    # :field s_end: 1-based number of subject sequence, in which alignment ends;
    # :type s_end: int;
    # :method __init__: accepts: (q_align, s_align, q_start, q_end, s_start, s_end);
    #     After calling the __init__ methof all these arguments will
    #     become corresponding fields of class's instance
    # :method get_align_len: returns length of alignment length of 'int';
    # :method get_pident: returns identity of the alignment of 'float' [0, 1];

    def __init__(self, q_align, s_align, q_start, q_end, s_start, s_end, qlen,
        sacc, sstrand, score, length=None, gaps=None, evalue=None):
        # :param q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
        # :type q_align: str;
        # :param s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
        # :type s_align: str;
        # :param q_start: 1-based number of query sequence, in which alignment starts;
        # :type q_start: int;
        # :param q_end: 1-based number of query sequence, in which alignment ends;
        # :type q_end: int;
        # :param s_start: 1-based number of subject sequence, in which alignment starts;
        # :type s_start: int;
        # :param s_end: 1-based number of subject sequence, in which alignment ends;
        # :type s_end: int;


        self.q_align = q_align
        self.s_align = s_align

        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end = s_end

        self.qlen = qlen

        if not length is None:
            self.length = length
        else:
            self.length = len(q_align)
        # end if

        if not gaps is None:
            self.gaps = gaps
        else:
            self.gaps = round(q_align.count('-') / self.length, 2) 
        # end if

        self.sacc = sacc
        self.sstrand = sstrand

        self.score= score
        self.evalue = evalue
    # end def __init__

    def __repr__(self):
        return "\nQS={}; QE={}; SS={}; SE={}".format( self.q_start, self.q_end,
                                                            self.s_start, self.s_end )
    # end def __repr__
# end class AlignResult


def SW_get_new_score(up, left, diag,
                  matched, gap_penalty, match, mismatch):
    # Function is dedicated to fill an element of alignment matrix during
    # forward step of Smith-Waterman algorithm.
    # :param up: value of the upper element of matrix;
    # :type up: int;
    # :param left: value of the left element of matrix;
    # :type left: int;
    # :param diag:value of the diagonal element of matrix;
    # :type diagonal: int;
    # :param matched: logic value: is True if nucleotides are equal and False otherwise;
    # :type matched: bool;
    # :param gap_penalty: penalty for gap;
    # :type gap_penalty: int;
    # :param match: reward for matching;
    # :type match: int;
    # :param mismatch: penalty for mismatch;
    # :type mismatch: int;

    add_to_diag = match if matched else mismatch
    return max(0, diag+add_to_diag, left-gap_penalty, up-gap_penalty)
# end def SW_get_new_score


# Constants for traceback.
END, DIAG, UP, LEFT = range(4)


def SW_next_move(SWsm, i, j, match, mismatch, gap_penalty, matched):
    # Function is dedicated to perform a move during traceback of Smith-Waterman algorithm.
    # :param SWsm: Smith-Waterman scoring matrix;
    # :type SWsm: list<list<int>>;   meant to be changed to list<array.array<unsigned_int>>
    # :param i: row index of scoring matrix;
    # :type i: int;
    # :param j: column index of scoring matrix;
    # :type j: int;
    # :param gap_penalty: penalty for gap;
    # :type gap_penalty: int;
    # :param matched: logic value: is True if nucleotides are equal and False otherwise;
    # :type matched: bool;
    # Returns 1 if move is diagonal, 2 if move is upper, 3 if move is left.
    # Raises SWInvalidMove when Smith-Waterman algorithm can't perform the next move.

    diag = SWsm[i - 1][j - 1]
    up   = SWsm[i - 1][  j  ]
    left = SWsm[  i  ][j - 1]

    if (diag + (mismatch, match)[matched]) == SWsm[i][j]:
        return DIAG if diag!=0 else END
    
    elif (up - gap_penalty) == SWsm[i][j]:
        return UP if up!=0 else END
    
    elif (left - gap_penalty) == SWsm[i][j]:
        return LEFT if up!=0 else END
    
    else:
        # Execution should not reach here.
        raise SWInvalidMove("""Smith-Waterman algorithm: invalid move during traceback:
            diag={}; up={}; left={}; matched={}; SWsm[i][j]={}""".format(diag, up, left, matched, SWsm[i][j]))
    # end if
# end def SW_next_move


def SW_align(jseq, iseq, sacc, gap_penalty=5, match=1, mismatch=-1):
    # Function performs Smith-Waterman aligning algorithm.
    # Reward and all penalties are int for the sake of acceleration of this slow snake.
    # :param jseq: 'query' sequence, which is loated to the top of scoring matrix, across the j index varying;
    # :type jseq: str;
    # :param iseq: 'subject' sequence, which is loated to the left of scoring matrix, across the i index varying;
    # :type iseq: str;
    # :param gap_penalty: penalty for gap;
    # :type gap_penalty: int;
    # :param match: reward for matching;
    # :type match: int;
    # :param mismatch: penalty for mismatch;
    # :type mismatch: int;
    # Returns instance of AlignResult performing results of alignment.

    SWsm = [array('I', [0 for j in range(len(jseq)+1)]) for i in range(len(iseq)+1)]
    # SWsm = [[0 for j in range(len(jseq)+1)] for i in range(len(iseq)+1)]
    score = partial(SW_get_new_score, gap_penalty=gap_penalty, match=match, mismatch=mismatch)
    max_val, max_pos = 0, None

    for i in range(1, len(iseq)+1):
        for j in range(1, len(jseq)+1):
            SWsm[i][j] = score(up=SWsm[i-1][j], left=SWsm[i][j-1], diag=SWsm[i-1][j-1], matched=iseq[i-1] == jseq[j-1])

            if SWsm[i][j] > max_val:
                max_val = SWsm[i][j]
                max_pos = (i, j)
            # end if
        # end for
    # end for

    # It happens (no alignment)
    if max_pos is None:
        return None
    # end if

    aligned_jseq = ""
    aligned_iseq = ""
    i, j = max_pos
    move = SW_next_move(SWsm, i, j, match, mismatch,
                        gap_penalty, iseq[i-1] == jseq[j-1])

    while move != END:

        if move == DIAG:
            aligned_jseq += jseq[j-1]
            aligned_iseq += iseq[i-1]
            j -= 1
            i -= 1

        elif move == UP:
            aligned_jseq += '-'
            aligned_iseq += iseq[i-1]
            i -= 1

        elif move == LEFT:
            aligned_jseq += jseq[j-1]
            aligned_iseq += '-'
            j -= 1
        # end if

        move = SW_next_move(SWsm, i, j, match, mismatch,
                            gap_penalty, iseq[i-1] == jseq[j-1])
    # end while

    if jseq[j-1] == iseq[i-1]:
        aligned_jseq += jseq[j-1]
        aligned_iseq += iseq[i-1]
    # end if

    return AlignResult(aligned_jseq[::-1], aligned_iseq[::-1], j, max_pos[1], i, max_pos[0], len(jseq),
        sacc, True, max_val)
# end def SW_align
