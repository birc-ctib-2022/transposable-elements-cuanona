[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9399879&assignment_repo_type=AssignmentRepo)
# Project 5: Simulating transposable elements

In the last project, we imagine that someone has hired us to help out with simulating a genome containing [transposable elements]. (I know people who have such strange interests, so it is not beyond the realm of possibilities).

We won’t do anything complicated, this is just an exercise after all, but we will want to simulate TEs as stretches of DNA that can copy themselves elsewhere in the genome.

Our employer already has most of the simulator up and running. She has a program that randomly picks operations to do—insert a TE ab initio, copy a TE, or disable one with a mutation—but she needs us to program a representation of a genome to track where the TEs are.

There are multiple ways to do this, but you should implement at least two: one-based Python lists, where each nucleotide is represented by one entry in a list, and one based on linked lists, where each nucleotide is represented by a link. If you feel ambitious, you can try others (for example keeping track of ranges of a genome with the same annotation so you don’t need to explicitly represent each nucleotide).

## Genome interface

A genome should be represented as a class that implements the following methods:

```python
class Genome(ABC):
    """Representation of a circular genome."""

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
        the copyleft of index 0 or the right of the largest index, it should
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
        linear, but imagine that the last character is immediately followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        ...  # not implemented yet

```

The `ABC` and `@abstractmethod` just means that this class is not something you can use by itself, but that another class must implement the details. In `src/genome.py` you will find templates for a Python list and a linked list implementation (without the actual implementation, because you have to implement them).

You are free to implement the genome classes however you want, and use whatever auxiliary data structures you desire, as long as one uses a Python list with an element for each nucleotide and the other a linked list with a link for each nucleotide. If you want to implement a third (or fourth or fifth...) version, you are very welcome to do so as well.

## Complexity

When you have implemented the two (or more) classes, describe the complexity of each operation as a function of the genome size (at the time of the operation), and the size of the TE involved (and when copying, the offset you are copying). Put the description here:

## List implementation

During this exercise, we will use $n$ as the length of the genome, $m$ as the length of the transposable element, and $k$ as the number of "blocks". 

First, let's start with the List implementation. 

### Init instance

```python
def __init__(self, n: int):
    self.genome = n*['-'] # O(n)
    self.active_identifiers = {} # O(1)
    self.te_counter = 0# O(1)
```

Since appending to a list has a complexity of $O(1)$, creating the genome costs $O(n)$, whereas $n$ is the length of the initial genome. 

### Insert a transposable element to the list

```python
def insert_te(self, pos: int, length: int) -> int:
    pos = pos % len(self) # O(1)
    active_identifiers = list(self.active_identifiers.items()) # O(k)
    for identifier,active_te  in active_identifiers: # O(k)
        is_in_range(pos, active_te): # O(1)
            start, end = active_te # O(1)
            self.genome[start:end] = (end - start)*["x"] # O(n + m)
            self.active_identifiers.pop(identifier) # O(1)
        if active_te.start > pos: O(1)
            self.active_identifiers[identifier] = Interval(
                active_te.start + length, active_te.end + length
            ) # O(1)
    self.genome[pos:pos] = length*['A'] # O(n + m)
    self.te_counter += 1 # O(1)
    self.active_identifiers[self.te_counter] = Interval(pos, pos + length) # O(1)
    return self.te_counter
```

As you can see in the previous block, the final complexity is $O(k(n+m))$. However, in the worst case, all transposable elements are of length one and $k = n$, then we have $O(n^2\cdot m)$

### Copy a transposable element

```python
def copy_te(self, te: int, offset: int) -> int | None:
    original: Interval = self.active_identifiers.get(te) # O(1)
    if original: O(1)
        pos = (original.start + offset) # O(1)
        length = original.end - original.start # O(1)
        return self.insert_te(pos, length) # inserting transposable element
```

As we are using a dict (hash table) for storing the exact position of each active transposon, this operation is, after calculating the new position which can be done in $O(1)$, just inserting a transposable element. Then, in the worst case its complexity is the same as inserting an element and in the best case is $O(1)$. 

### Disable a transposable element

```python
def disable_te(self, te: int) -> None:
    original: Interval = self.active_identifiers.pop(te) # O(1)
    if original: O(1)
        self.genome[original.start:original.end] = (original.end - original.start)*["x"] # O(n+m)
```

Again, we can find the right indexes in $O(1)$. For disabling the te, we have to set a slice which, as [TimeComplexity Wiki](https://wiki.python.org/moin/TimeComplexity) says, its amortized worst case is $O(n, m)$.

## Get active in te

```python
def active_tes(self) -> list[int]:
    return list(self.active_identifiers.keys())
```

As we save all active te and update those entries in a dict in the other operations, getting that list is just $O(1)$. Or, if we copy it, just $O(k)$. 

## Get length of the genome

```python
def __len__(self) -> int:
    return len(self.genome)
```
As we store the length of the list with the list itself, this operation is $O(1)$. 

### Represent genome

```python
def __str__(self) -> str:
    return "".join(self.genome)
```

This operation is $O(n)$ because we are appending nucleotides to the str. 

## Linked list implementation with blocks

### Init instance

```python
def __init__(self, n: int):
    self.head = Link(None, None, None)  # type: ignore
    self.head.prev = self.head
    self.head.next = self.head
    insert_after(self.head.prev, Feature(self.empty_te, n))
```

All operations are $O(1)$. Notice that it doesn't depend on $n$ because we are considering blocks instead of nucleotides. 

### Insert a transposable element to the list

```python
def insert_te(self, pos: int, length: int) -> int:
    for feature, end in self.into_iter_with_pos(): # O(k)
        if end > pos: # O(1)
            self.insert_into(length, feature, feature.val.length -end + pos, end - pos) # O(1)
            return self.counter_te
    return -1
```

As you can see in the previous block, the final complexity is $O(k)$. However, in the worst case, all transposable elements are of length one and $k = n$, then we have $O(n)$. Notice that, as we are using blocks instead of nucleotides, the size of the transposon doesn't affect it. We iterate over, at most, $k$ blocks and insert the new node in constant time once we find the right place. 

### Copy a transposable element


```python
# This is simplified code
def copy_te(self, te: int, offset: int) -> int | None:
    original = feature= self.active_identifier[te] # O(1)
    if original.val.feature != self.active_te: #O(1)
        return None
    # Calculate offset -> O(1)
    # Iterate over the linked list until finding the node where to insert O(k) worst case
    while offset > 0:
        feature = getattr(feature, direction)
        if feature is self.head:
            feature = getattr(feature, direction)
            offset -= feature.val.length
    # Insert te
    return self.counter_te
```

For this implementation, we have to find the node to insert our new transposable element that has a complexity of $O(k)$. In the worst case, where all elements in our linked list have length one, that would be $O(n)$. After finding the node, we have to insert it. Then, the final complexity is $O(2k) = O(k)$.

### Disable a transposable element

```python
def disable_te(self, te: int) -> None:
    feature = self.active_identifier.pop(te) # O(1)
    if feature is not None: O(1)
        feature.val = Feature(self.inactive_te, feature.val.length) # O(1)
```

We can find the pointer to the element we want to modify in $O(1)$ and overwrite the tuple in $O(1)$ too. Then, its complexity is $O(1)$.

## Get active in te

```python
def active_tes(self) -> list[int]:
    return list(self.active_identifiers.keys())
```

As we save all active te and update those entries in a dict in the other operations, getting that list is just $O(1)$. Or, if we copy it, just $O(k)$. 

## Get the length of the genome

```python
def __len__(self) -> int:
    end = 0 #O(1)
    for _, end in self.into_iter_with_pos(): # O(k)
        pass
    return end
```
As we have to iterate over all elements until we have reached the first dummy element, this operation is $O(k)$. 

### Represent genome

```python
def __str__(self) -> str:
    mapping = "-Ax"
    return "".join(
        (el.length )*mapping[el.feature] for el in self
    )
```

This operation is $O(n)$ because we are appending all the nucleotides. 



In `src/simulate.py` you will find a program that can run simulations and tell you the actual time it takes to simulate with different implementations. You can use it to test your analysis. You can modify the parameters of the simulator if you want to explore how they affect the running time.
