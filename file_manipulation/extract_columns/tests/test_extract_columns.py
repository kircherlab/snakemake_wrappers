from click.testing import CliRunner
import filecmp

from extract_columns import cli


def test1_extract_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A'])
    assert result.exit_code == 2


def test2_extract_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 2


def test3_extract_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest3.tsv")


def test4_extract_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'B',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest4.tsv")
