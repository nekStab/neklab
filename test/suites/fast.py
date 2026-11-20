"""Level 2 -- fast tests (to be written).

Same cases as the full level, but each solver starts next to the answer of
a full run (converged base flow or leading mode + a little seeded noise), and
the test checks it gets back to the reference answer in few iterations.

Each fast test will be a complete case folder test/fast/<case>/, exactly like
the unit tests: its .usr loads the converged field, adds the noise, runs the
solver and prints NEKLAB_RESULT lines (e.g. newton_iterations,
eigs_matvecs) that the test compares with bounds and references.
"""
