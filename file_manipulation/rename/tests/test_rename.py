from click.testing import CliRunner
import filecmp

from rename import cli


def test1_rename():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 2


def test2_rename():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--row', 'A',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 2


def test3_rename():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', 'B',
                                 '--column', 'B', 'A',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest3.tsv")


def test4_rename():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', 'G',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest4.tsv")
