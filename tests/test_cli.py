"""Unit tests for GenoHMM CLI."""

from geno_hmm.cli import main


def test_cli_demo() -> None:
    code = main(["demo", "--scenario", "digeorge", "--noise", "0.8"])
    assert code == 0


def test_cli_benchmark() -> None:
    code = main(["benchmark", "--num-runs", "2"])
    assert code == 0
