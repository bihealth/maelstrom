from maelstrom import snakemake as ms


def test_which():
    assert ms.which("python") is not None
    assert ms.which("python").endswith("python")
    assert ms.which("pythonX") is None
