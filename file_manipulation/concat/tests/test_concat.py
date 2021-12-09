from click.testing import CliRunner
import filecmp

from concat import cli


def test1_concat():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv'])
    assert result.exit_code == 2


def test2_concat():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--column', '1',
                                 '--output', "tests/output.tsv"])
    assert result.exit_code == 2


def test3_concat():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--input', 'tests/input1.tsv',
                                 '--column', 'x', '1',
                                 '--column', 'x', '2',
                                 '--output', "tests/output3.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output3.tsv", "tests/outTest3.tsv")


def test4_concat():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input1.tsv',
                                 '--input', 'tests/input2.tsv',
                                 '--input', 'tests/input3.tsv',
                                 '--index', 'metric',
                                 '--output', "tests/output4.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output4.tsv", "tests/outTest4.tsv")
