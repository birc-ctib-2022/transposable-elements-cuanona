from collections import namedtuple

from genome import Genome


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
        index = self.find_where_to_insert(pos)
        self.disable_if_active(index)
        self.insert_into(self.active_te, index, pos, pos + length)
        self.update_te_counter()
        self.identifiers_active[self.counter_te] = index +1
        self.update_identifiers_from(pos)
        return self.counter_te

    def update_te_counter(self):
        self.counter_te += 1

    def disable_if_active(self, index):
        for identifier, i in self.identifiers_active.items():
            if i == index:
                self.disable_te(identifier)
                break

    def update_identifiers_from(self, pos):
        for k,value in self.identifiers_active.items():
            if value > pos:
                self.identifiers_active[k] += 3
    def find_where_to_insert(self, start: int):
        for index, stack in enumerate(self.genome):
            if stack.start <= start <= stack.end:
                return index
    def insert_into(self, new_feature: str, index: int, start: int, end: int):
        old: Feature = self.genome[index]
        split_size = old.end - start
        if split_size < 0:
            raise IndexError
        self.genome[index] = Feature(new_feature, start, end)
        self.genome.insert(
                    index,
                    Feature(old.feature, old.start, start)
                    )
        if index + 2 > len(self.genome):
            self.genome.append(Feature(old.feature, end, end + split_size))
        else:
            self.genome.insert(
                    index + 2,
                    Feature(old.feature, end, end + split_size)
                    )
        self.update_genome_by_add_diff(end - start, index+3)

    def update_genome_by_add_diff(self, length, index):

        for i in range(index, len(self.genome)):
            old:Feature = self.genome[i]
            self.genome[i] = Feature(
                old.feature,
                old.start + length,
                old.end + length
                )    
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
        original_index = self.identifiers_active.get(te)
        if original_index is None:
            return None
        original = self.genome[original_index]     
        new_position = self.calculate_pos_from_offset(offset, original)
        return self.insert_te(new_position, original.end - original.start)

    def calculate_pos_from_offset(self, offset, original):
        new_position = original.start + offset
        match new_position:
            case _ if new_position < 0:
                return len(self) + new_position
            case _ if new_position > len(self):
                return new_position -len(self)
        return new_position

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        original_index = self.identifiers_active.get(te)
        if original_index:
            self.identifiers_active.pop(te)
            original = self.genome[original_index]  
            self.genome[original_index] = Feature(
                self.inactive_te,
                original.start,
                original.end
            )

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.identifiers_active.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        return self.genome[-1].end

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
