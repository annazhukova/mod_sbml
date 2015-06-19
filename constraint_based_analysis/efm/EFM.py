import math
import sys
import gmpy

__author__ = 'anna'


def get_int_size():
    """
    Calculates the maximal number of bits in an int:
    math.log(sys.maxint) / math.log(2).

    :return: int, the maximal number of bits in an int.
    """
    return math.log(sys.maxint) / math.log(2)


def get_binary_efm_len(binary_efm):
    return sum(gmpy.popcount(it) for it in binary_efm)


class EFM(object):
    def __init__(self, r_ids, rev_r_ids=None, int_size=None, r_id2coeff=None, binary_efm=(), coefficients=None):
        """
        Converts an EFM representation {r_id: coefficient} into an inner binary representation
        (binary_representation, non_zero_coefficients).

        A binary representation of an EFM is a list of integers whose binary representations
        correspond to the reactions that are active in the EFM: if the reaction is active,
        the corresponding bit is set to 1.
        If the total number of reactions in the model is larger that the number of bits in an int,
        several ints are used.

        Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
        a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
        that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
        The non_zero_coefficients for this EFM are [3, -2, 1, 1].

        :param r_id2coeff: a EFM represented as a dictionary {r_id: coefficient}.

        :param r_ids: ordered collection of reaction ids (strings).

        :param rev_r_ids: (optional) set of ids of reversible reactions (strings).
        If not present all reactions are assumed to be irreversible.

        :param int_size: (optional) the maximal number of bits in an int,
        can be calculated as math.log(sys.maxint) / math.log(2).
        """
        if (not r_id2coeff and not binary_efm) or not r_ids:
            raise AttributeError('Either r_id2coeff or binary_efm, and r_ids should be specified')
        else:
            if not rev_r_ids:
                rev_r_ids = set()
            self.r_ids = r_ids
            self.rev_r_ids = rev_r_ids
            self.int_size = int_size if int_size else get_int_size()
            self.binary_efm = binary_efm
            self.coefficients = coefficients
            if r_id2coeff:
                self.__from_r_id2coeff(r_id2coeff)

    def __str__(self, binary=False):
        return self.to_string()

    def to_string(self, binary=False):
        r_id2coefficient = self.to_r_id2coeff()
        keys = sorted(r_id2coefficient.iterkeys())
        if binary or not self.coefficients:
            return '\t'.join('%s%s' % ('-' if r_id2coefficient[r_id] < 0 else '', r_id) for r_id in keys)
        return '\t'.join('%g %s' % (r_id2coefficient[r_id], r_id) for r_id in keys)

    def __from_r_id2coeff(self, r_id2coeff):
        """
        Converts an EFM representation {r_id: coefficient} into a tuple (binary_representation, non_zero_coefficients).

        A binary representation of an EFM is a list of integers whose binary representations
        correspond to the reactions that are active in the EFM: if the reaction is active,
        the corresponding bit is set to 1.
        If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

        Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
        a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
        that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
        The non_zero_coefficients for this EFM are [3, -2, 1, 1].

        :param r_id2coeff: a EFM represented as a dictionary {r_id: coefficient}.
        """
        bit_efm = []
        bit_efm_cur = 0
        i = 1
        coefficients = []

        def advance(i, bit_efm_cur):
            if math.log(i) == self.int_size:
                bit_efm.append(bit_efm_cur)
                bit_efm_cur = 0
                i = 1
            else:
                i <<= 1
            return bit_efm_cur, i

        for r_id in self.r_ids:
            if r_id in r_id2coeff:
                coeff = r_id2coeff[r_id]
            else:
                coeff = 0
            if coeff > 0:
                bit_efm_cur |= i
                coefficients.append(coeff)
            bit_efm_cur, i = advance(i, bit_efm_cur)
            if r_id in self.rev_r_ids:
                if coeff < 0:
                    bit_efm_cur |= i
                    coefficients.append(coeff)
                bit_efm_cur, i = advance(i, bit_efm_cur)
        bit_efm.append(bit_efm_cur)

        self.binary_efm = tuple(bit_efm)
        self.coefficients = tuple(coefficients)

    def to_r_id2coeff(self, binary=False):
        """
        Returns a representation of this EFM as a dictionary
        that maps ids of active reactions to their coefficients: {r_id: coefficient}.
        If binary is set to True, the coefficient values are 1 (for any active reaction in the standard direction)
        or -1 (for reactions that are active in the reversed direction).

        :param binary: boolean, if is set to True, the coefficient values in the result will be
        1 (for any active reaction in the standard direction)
        or -1 (for reactions that are active in the reversed direction).
        Otherwise, any float coefficient values will be possible.

        :return: dict, {r_id: coefficient}.
        """
        converted_efm = {}
        binary_efm_iterable = iter(self.binary_efm)
        cur_efm_part = next(binary_efm_iterable)
        coeff_iterable = iter(self.coefficients) if self.coefficients else None

        def process(r_id, cur_efm_part, i, reversed=False):
            if cur_efm_part & i:
                coeff = next(coeff_iterable) if self.coefficients else (-1 if reversed else 1)
                if binary:
                    converted_efm[r_id] = 1 if coeff > 0 else -1
                else:
                    converted_efm[r_id] = coeff
            if math.log(i) == self.int_size:
                cur_efm_part = next(binary_efm_iterable)
                i = 1
            else:
                i <<= 1
            return i, cur_efm_part

        i = 1
        for r_id in self.r_ids:
            i, cur_efm_part = process(r_id, cur_efm_part, i)
            if r_id in self.rev_r_ids:
                i, cur_efm_part = process(r_id, cur_efm_part, i, True)
        return converted_efm

    def intersection(self, other):
        if not other or not isinstance(other, EFM):
            raise AttributeError('Other should be of type EFM')
        return EFM(binary_efm=tuple((p1 & p2 for (p1, p2) in zip(self.binary_efm, other.binary_efm))),
                   r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, int_size=self.int_size)

    def __len__(self):
        """
        Returns the length of this EFM, i.e. the number of active reactions.
        """
        return len(self.coefficients) if self.coefficients else get_binary_efm_len(self.binary_efm)

    def __eq__(self, other):
        if not other or not isinstance(other, EFM):
            return False
        return (self.binary_efm, self.coefficients) == (other.binary_efm, other.coefficients)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.binary_efm, self.coefficients))



