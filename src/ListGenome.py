from collections import namedtuple

from genome import Genome


"""A circular genome for simulating transposable elements."""

Range = namedtuple("Range", ["start", "end"])

def is_in_range(pos: int, range: Range):
    return range.start <= pos < range.end

class ListGenome(Genome):
    """Representation of a circular enome."""
    genome: list[str]
    active_identifiers: dict[int:Range]
    def __init__(self, n: int):
        """Create a genome of size n."""
        self.genome = n*['-']
        self.active_identifiers = {}
        self.te_counter = 0

    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        pos = pos % len(self)
        active_identifiers = list(self.active_identifiers.items())
        for identifier,active_te  in active_identifiers:
            if is_in_range(pos, active_te):
                start, end = active_te
                self.genome[start:end] = (end - start)*["x"]
                self.active_identifiers.pop(identifier)
                pass
            if active_te.start > pos:
                self.active_identifiers[identifier] = Range(
                    active_te.start + length, active_te.end + length
                )
        self.genome = self.genome[:pos] + length*['A'] + self.genome[pos:]
        self.te_counter += 1
        self.active_identifiers[self.te_counter] = Range(pos, pos + length)
        return self.te_counter

    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        original: Range = self.active_identifiers.get(te)
        if original:
            pos = (original.start + offset)
            length = original.end - original.start
            return self.insert_te(pos, length)

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        original: Range = self.active_identifiers.pop(te)
        if original:
            self.genome[original.start:original.end] = (original.end - original.start)*["x"]

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.active_identifiers.keys())

    def __len__(self) -> int:
        """Get the current length of the genome."""
        return len(self.genome)

    def __str__(self) -> str:
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        return "".join(self.genome)


