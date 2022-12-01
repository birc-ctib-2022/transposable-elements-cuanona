
from __future__ import annotations
from typing import (
    Generic, TypeVar, Iterable,
    Callable, Protocol
)

from collections import namedtuple

from genome import Genome


Feature = namedtuple("Feature", ["feature", "length"])

T = TypeVar('T')

class Link(Generic[T]):
    """Doubly linked link."""

    val: T
    prev: Link[T]
    next: Link[T]

    def __init__(self, val: T, p: Link[T], n: Link[T]):
        """Create a new link and link up prev and next."""
        self.val = val
        self.prev = p
        self.next = n

def insert_after(link: Link[T], val: T) -> None:
    """Add a new link containing avl after link."""
    new_link = Link(val, link, link.next)
    new_link.prev.next = new_link
    new_link.next.prev = new_link
def insert_before(link: Link[T], val: T) -> None:
    """Add a new link containing avl after link."""
    new_link = Link(val, link, link.prev)
    new_link.prev.next = new_link
    new_link.next.prev = new_link

class LinkedListGenome(Genome):
    """
    Representation of a genome.

    Implements the Genome interface using linked lists.
    """
    head: Link[Feature]  # Dummy head link
    empty_te, active_te, inactive_te = 0, 1, 2
    active_identifier = {}
    counter_te = 0

    def __init__(self, n: int):
        """Create a new genome with length n."""
        self.head = Link(None, None, None)  # type: ignore
        self.head.prev = self.head
        self.head.next = self.head
        insert_after(self.head.prev, Feature(self.empty_te, n))
    
    def __iter__(self):
        feature = self.head.next
        while feature is not self.head:
            yield feature.val
            feature = feature.next
    def into_iter_with_pos(self):
        feature = self.head.next
        acc = 0
        while feature is not self.head:
            acc += feature.val.length 
            yield (feature, acc)
            feature = feature.next

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

        for feature, end in self.into_iter_with_pos():
            if end > pos:
                self.insert_into(length, feature, feature.val.length -end + pos, end - pos)
                return self.counter_te
        return -1

    def insert_into(self, length, feature, split_size_1, split_size_2):
        if  feature.val.feature == self.active_te:
            self.disable_feature(feature)
        feature.val = Feature(feature.val.feature, split_size_1)
        insert_after(feature, Feature(self.active_te, length))
        insert_after(feature.next, Feature(feature.val.feature, split_size_2))
        self.counter_te += 1
        self.active_identifier[self.counter_te] = feature.next

    def disable_feature(self, feature):
        feature.val = Feature(self.inactive_te, feature.val.length)
        for identifier, disable in self.active_identifier.items():
            if disable is feature:
                self.active_identifier.pop(identifier)
                break

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
        original = feature= self.active_identifier[te]
        if original.val.feature != self.active_te:
            return None
        match offset > 0:
            case True:
                direction = "next"
                offset = offset - original.val.length
            case False:
                direction = "prev"
                offset *= -1
        while offset > 0:
            feature = getattr(feature, direction)
            if feature is self.head:
                feature = getattr(feature, direction)
            offset -= feature.val.length
        match direction:
            case "next":
                self.insert_into(
                    original.val.length,
                    feature,
                    feature.val.length + offset,
                    abs(offset)
                    )
            case "prev":
                self.insert_into(
                    original.val.length, feature,
                    abs(offset),
                    feature.val.length + offset
                    )
        return self.counter_te

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        feature = self.active_identifier.pop(te)
        if feature is not None:
            feature.val = Feature(self.inactive_te, feature.val.length)

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.active_identifier.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        end = 0
        for _, end in self.into_iter_with_pos():
            pass
        return end


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
            (el.length )*mapping[el.feature] for el in self
            )