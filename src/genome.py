"""A circular genome for simulating transposable elements."""

from abc import (
    # A tag that says that we can't use this class except by specialising it
    ABC,
    # A tag that says that this method must be implemented by a child class
    abstractmethod
)
from collections import namedtuple
Feature = namedtuple("Feature", ["feature", "start", "end"])

class Genome(ABC):
    """Representation of a circular enome."""

    def __init__(self, n: int):
        """Create a genome of size n."""

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """

    @abstractmethod
    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""

    @abstractmethod
    def __len__(self) -> int:
        """Get the current length of the genome."""

    @abstractmethod
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

Feature = namedtuple("Feature", ["feature", "start", "end"])

class ListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using Python's built-in lists
    """

    genome: list[Feature]
    identifiers_active: dict[int, int]

    empty_te, active_te, inactive_te = 0, 1, 2
    

    def __init__(self, n: int):
        """Create a new genome with length n."""
        self.genome = [Feature(self.empty_te, 0, n)]
        self.identifiers_active = {}
        self.counter_te = 0
        self.identifiers_inactive = list()

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
        new_te = Feature(self.active_te, pos, pos + length)
        for index, old_te in enumerate(self.genome):
            if old_te.end > pos:
                diff = old_te.end - pos
                return self.split_and_insert(new_te, old_te, index, diff)
                
    def split_and_insert(self, new_te, old_te, index, diff):
        """Split and insert a new te"""
        if old_te.feature == self.active_te:
            old_te = Feature(
                        self.inactive_te,
                        old_te.start,
                        old_te.end
                    )
        self.genome[index] = new_te
        self.genome.insert(
                    index,
                    Feature(old_te.feature, old_te.start, old_te.end - diff)
                    )
        self.genome.insert(
                    index + 2,
                    Feature(old_te.feature, new_te.end, new_te.end + diff)
                    )
        self.counter_te += 1
        new_identifiers = {}
        for identifier, dict_pos in self.identifiers_active.items():
            match dict_pos:
                case dict_pos  if dict_pos < index:
                    new_identifiers[identifier] = dict_pos
                case dict_pos if dict_pos > index:
                     new_identifiers[identifier] = dict_pos + 2
        new_identifiers[self.counter_te] = index + 1
        self.identifiers_active = new_identifiers
        return self.counter_te
        
    
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
        original_index = self.identifiers_active[te]
        original = self.genome[original_index]
        original_position = original.start
        new_position = original_position + offset
        return self.insert_te(new_position, original.end - original.start)
        

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # FIXME

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.identifiers_active.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        ...  # FIXME
        return 0

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
        mapping = "-Ax"
        return "".join(
            (el.end - el.start)*mapping[el.feature] for el in self.genome
            )


class LinkedListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using linked lists.
    """

    def __init__(self, n: int):
        """Create a new genome with length n."""
        

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
        ...  # FIXME
        return -1

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
        ...  # FIXME

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # FIXME

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        # FIXME
        return []

    def __len__(self) -> int:
        """Current length of the genome."""
        # FIXME
        return 0

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
        return "FIXME"
