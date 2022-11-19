# Project 5: Simulating transposable elements

In the last project, we imagine that someone has hired us to help out with simulating a genome containing [transposable elements]. (I know people who has such strange interests, so it is not beyond the realm of possibilities).

We won’t do anything complicated, this is just an exercise after all, but we will want to simulate TEs as stretches of DNA that can copy themselves elsewhere in the genome.

Our employer already has most of the simulator up and running. She has a program that randomly picks operations to do—insert a TE ab initio, copy a TE, or disable one with a mutation—but she needs us to program a representation of a genome to track where the TEs are.

There are multiple ways to do this, but you should implement at least two: one based Python lists, where each nucleotide is represented by one entry in a list, and one based on linked lists, where each nucleotide is represented by a link. If you feel ambitious, you can try others (for example keeping track of ranges of a genome with the same annotation so you don’t need to explicitly represent each nucleotide).

## Genome interface

A genome should be represented as a class that implements the following methods:

```python
class Genome(ABC):
    """Representation of a circular enome."""

    def __init__(self, n: int):
        """Create a genome of size n."""
        ...  # not implemented yet

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
        ...  # not implemented yet

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
        ...  # not implemented yet

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # not implemented yet

    @abstractmethod
    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        ...  # not implemented yet

    @abstractmethod
    def __len__(self) -> int:
        """Get the current length of the genome."""
        ...  # not implemented yet

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
        ...  # not implemented yet

```

The `ABC` and `@abstractmethod` just means that this class is not something you can use by itself, but that another class must implement the details. In `src/genome.py` you will find templates for a Python list tand a linked list implementation (without the actual implementation, because you have to implement them).

You are free to implement the genome classes however you want, and using whateer auxilary data structures you desire, as long as one uses a Python list with an element for each nucleotide and the other a linked list with a link for each nucleotide. If you want to implement a third (or fourth or fifth...) version, you are very welcome to do so as well.

## Complexity

When you have implemented the two (or more) classes, describe the complexity of each operation as a function of the genome size (at the time of the operation), and the size of the TE involved (and when copying, the offset you are copying). Put the description here:

**FIXME: OPERATION COMPLEXITY**

In `src/simulate.py` you will find a program that can run simulations and tell you actual time it takes to simulate with different implementations. You can use it to test your analysis. You can modify the parameters to the simulator if you want to explore how they affect the running time.
