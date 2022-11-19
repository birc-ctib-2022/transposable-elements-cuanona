"""Testing genome implementations."""


# This directory will be checked with pytest. It will examine
# all files that start with test_*.py and run all functions with
# names that start with test_

from genome import (
    Genome,
    ListGenome,
    LinkedListGenome
)
from typing import Type


def run_genome_test(genome_class: Type[Genome]) -> None:
    """Test a Genome implementation."""
    genome = genome_class(20)
    assert str(genome) == "--------------------"
    assert genome.active_tes() == []

    assert 1 == genome.insert_te(5, 10)   # Insert te 1
    assert str(genome) == "-----AAAAAAAAAA---------------"
    assert genome.active_tes() == [1]

    assert 2 == genome.insert_te(10, 10)  # Disable 1 but make 2 active
    assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx---------------"
    assert genome.active_tes() == [2]

    # Make TE 3 20 to the right of the start of 2
    assert 3 == genome.copy_te(2, 20)
    assert str(genome) == "-----xxxxxAAAAAAAAAAxxxxx-----AAAAAAAAAA----------"
    assert genome.active_tes() == [2, 3]

    # Make TE 4 15 to the leftt of the start of 2
    assert 4 == genome.copy_te(2, -15)
    assert str(genome) \
        == "-----xxxxxAAAAAAAAAAxxxxx-----AAAAAAAAAA-----AAAAAAAAAA-----"
    assert genome.active_tes() == [2, 3, 4]

    assert 5 == genome.insert_te(50, 10)
    assert str(genome) == \
        "-----xxxxxAAAAAAAAAAxxxxx-----" + \
        "AAAAAAAAAA-----xxxxxAAAAAAAAAAxxxxx-----"
    assert genome.active_tes() == [2, 3, 5]

    genome.disable_te(3)
    assert str(genome) == \
        "-----xxxxxAAAAAAAAAAxxxxx-----" \
        "xxxxxxxxxx-----xxxxxAAAAAAAAAAxxxxx-----"
    assert genome.active_tes() == [2, 5]


def test_list_genome() -> None:
    """Test that the Python list implementation works."""
    run_genome_test(ListGenome)


def test_linked_list_genome() -> None:
    """Test that the linked list implementation works."""
    run_genome_test(LinkedListGenome)
